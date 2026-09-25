import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, marching_cubes   
import pyvista
from . import fmmesh
from . import fmbeads
import pathlib

def tif_reader(input_file,color_idx):
	"""Intakes .tif files and returns strength of particular color across all voxels

	Parameters
	----------
	input_file : str
		String containing the path and filename format
		Example : './CytoD/Cell/Gel 2 CytoD%s.tif'
	color_idx : int
		The color to examine (0=red, 1=green, 2=blue)

	Returns 
	----------
	all_array : numpy.ndarray
		A NumPy array of shape (size_x,size_y,num_images) specifying the strength 
		of color_idx across the set of images

	"""

	parent_folder = str(pathlib.Path(input_file).parent)
	fnames = glob.glob(parent_folder + '/*.tif')
	num_images = len(fnames)
	sample_img = plt.imread(input_file%('0000'))
	size_x = sample_img.shape[0]
	size_y = sample_img.shape[1] 
	
	all_array = np.zeros((size_x,size_y,num_images))
	for kk in range(0,num_images):
		if kk < 10:
			num = '000%i'%(kk)
		elif kk < 100:
			num = '00%i'%(kk)
		else:
			num = '0%i'%(kk)
		
		fname =  input_file%(num)
		img = plt.imread(fname)
		all_array[:,:,kk] = img[:,:,color_idx]
	
	return all_array 

def get_cell_surface(input_file, dims, color_idx=0, cell_threshold=1.0):
	"""Creates an FMMesh object from image data

	Parameters
	----------
	input_file : str
		String containing the path and filename format
		Example : './CytoD/Cell/Gel 2 CytoD%s.tif'
	dims : np.array
		Total length of microscope imagery along the x, y, and z dimensions
	color_idx : int
		The color to examine (0=red, 1=green, 2=blue)
	cell_threshold : float
		Minimum voxel color intensity for consideration as part of the cell

	Returns 
	----------
	mesh : FMMesh
		An FMMesh object specifying the cell surface created from the image files

	"""

	X_DIM = dims[0]; Y_DIM = dims[1]; Z_DIM = dims[2]

	mesh = fmmesh.FMMesh()

	# import the image file and apply a gaussian filter 
	all_array = tif_reader(input_file,color_idx)
	all_array = ndimage.gaussian_filter(all_array,2)
	
	# threshold the image based on a set threshold
	bw = all_array > cell_threshold
	
	# find connected volumes
	label_img = label(bw, connectivity=bw.ndim)
	props = regionprops(label_img)
	centroids = np.zeros((len(props),3))
	areas = np.zeros((len(props)))
	for kk in range(len(props)):
		centroids[kk] = props[kk].centroid
		areas[kk] = props[kk].area
	
	# assume cell is the largest connected volume
	arg = np.argmax(areas)
	
	# save the cell volume
	vox_size = X_DIM / bw.shape[1] * Y_DIM / bw.shape[0] * Z_DIM / bw.shape[2]
	vol = areas[arg] * vox_size
	mesh.vol = vol
	
	# save the cell center
	cell_center = np.asarray([ centroids[arg,1] * X_DIM / bw.shape[1] , centroids[arg,0] * Y_DIM / bw.shape[0] , centroids[arg,2] * Z_DIM / bw.shape[2] ])
	mesh.center = cell_center
	
	# isolate the cell 
	bw_cell = np.zeros(bw.shape)
	for ii in range(0,bw.shape[0]):
		for jj in range(0,bw.shape[1]):
			for kk in range(0,bw.shape[2]):
				if label_img[ii,jj,kk] == int(arg+1):
					bw_cell[ii,jj,kk] = 1.0
	
	# flip the ii and jj dimensions to be compatible with the marching cubes algorithm
	bw_cell = np.swapaxes(bw_cell,0,1)
	
	# get the cell surface mesh from the marching cubes algorithm and the isolated cell image
	# https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.marching_cubes
	verts,faces, normals,_ = marching_cubes(bw_cell,spacing=(X_DIM/bw_cell.shape[0], Y_DIM/bw_cell.shape[1], Z_DIM/bw_cell.shape[2]))
	
	# save surface mesh info
	mesh.points = verts
	mesh.normals = normals
	mesh.faces = faces

	return mesh


##########################################################################################
# functions to pre-process bead images
#	OUTPUTS:
#		- x, y, z position of each bead based on the input images 
##########################################################################################
def signal_depth(all_array):
	"""Number of leading z slices that still carry bead signal.

	The whole-volume Otsu threshold is used as the test: a slice counts as
	carrying signal if its brightest voxel exceeds it. Unlike the per-slice
	threshold this one cannot adapt to an empty slice, so it is a meaningful
	statement about that slice's content.

	`all_array` is expected to be the *filtered* volume, so that the maxima
	and the threshold sit on the same scale as the per-slice Otsu that follows.

	The signal region is contiguous from the coverslip up in practice (checked
	over the 42 stacks of the ASMC dataset: for every one, the first slice
	below threshold is also the last slice above it). If a stray slice above
	the drop does clear the threshold, the returned depth extends to include
	it -- keeping a few noisy slices is the lesser error, since the caller
	still thresholds them per-slice.
	"""
	above = all_array.max(axis=(0,1)) > threshold_otsu(all_array)
	if not above.any():
		return 0
	return int(np.max(np.where(above)[0])) + 1


def get_bead_centers(input_file, dims, color_idx=1, threshold='per-slice',
		sigma=1):
	"""Creates a FMBeads object from image data

	Parameters
	----------
	input_file : str
		String containing the filename format
		Example : input_file='./CytoD/Beads/Gel 2 CytoD%s.tif'
	dims : np.array
		Total length of microscope imagery along the x, y, and z dimensions
	color_idx :
		The color to examine (0=red, 1=green, 2=blue)
	threshold : str or float
		Thresholding strategy; see `get_bead_centers_from_array`.
	sigma : float or sequence of 3 floats
		Gaussian smoothing before thresholding; see
		`get_bead_centers_from_array`.

	Returns
	----------
	beads : FMBeads
		An FMBeads object with bead positions corresponding to those calculated from imagery data

    """

	# import the image file
	all_array = tif_reader(input_file,color_idx)

	return get_bead_centers_from_array(all_array, dims, threshold=threshold,
		sigma=sigma)


def get_bead_centers_from_array(all_array, dims, threshold='per-slice',
		sigma=1):
	"""Creates a FMBeads object from an already-loaded intensity volume.

	Split out of `get_bead_centers` so that callers holding a volume from a
	reader other than `tif_reader` (e.g. an OME-TIFF z-stack) run the exact
	same segmentation. `get_bead_centers` is this function plus `tif_reader`.

	Parameters
	----------
	all_array : numpy.ndarray
		Intensity volume of shape (num_rows, num_cols, num_slices), i.e. the
		layout `tif_reader` returns: axis 0 is the image row (y), axis 1 the
		image column (x), axis 2 the z slice.
	dims : np.array
		Total length of microscope imagery along the x, y, and z dimensions
	threshold : str or float
		How the binary bead mask is thresholded, after the gaussian filter.

		'per-slice'
			An Otsu threshold recomputed for every z slice. The original
			behaviour, and the default.
		'global'
			One Otsu threshold over the whole volume.
		'per-slice-truncated'
			The global Otsu threshold is used only to find the depth at which
			the signal ends -- the first z slice whose maximum falls below it.
			That slice and everything above it are zeroed, and the per-slice
			Otsu is then run on what remains.
		float
			That absolute intensity, used directly.
	sigma : float or sequence of 3 floats
		Standard deviation, in voxels, of the Gaussian filter applied before
		thresholding. A sequence is in the array's own axis order (row, col,
		slice) = (y, x, z); 0 on an axis leaves that axis unsmoothed. The
		default 1 is the original isotropic filter.

		The filter keeps a bead's mask connected along z, where the
		point-spread function is several times longer than in x and y: on
		unsmoothed data a single slice falling below its Otsu threshold cuts
		the column in two and one bead is reported twice, stacked in z. But
		the same filter merges lateral neighbours closer than about 1.5 um
		(at 0.29 um voxels). On the iPSC stacks (0, 0, 1) keeps the z
		connectivity and resolves ~18% more beads than 1.

		Per-slice Otsu assumes every slice holds beads: Otsu maximizes the
		separation between two classes, so on a slice that is entirely
		background it splits the noise and returns a near-zero threshold,
		admitting anything nonzero. Where the bead signal dies partway up a
		stack -- which is the case for stacks acquired through an air
		objective into an aqueous sample -- the empty slices then contribute
		a large volume of thresholded noise, and `label(..., connectivity=3)`
		links it into a single object that swallows the real beads.

		A 'global' threshold removes that per-slice degree of freedom, but it
		is itself dragged down by all the empty slices it is averaged over, so
		in the bead-bearing slices it sits well below the per-slice value and
		merges neighbouring beads into blobs. 'per-slice-truncated' uses the
		global threshold for the one thing it is reliable for -- saying where
		the signal stops -- and leaves the thresholding inside the signal
		region to the per-slice Otsu that was calibrated for it.

	Returns
	----------
	beads : FMBeads
		An FMBeads object with bead positions corresponding to those calculated from imagery data

    """

	X_DIM = dims[0]; Y_DIM = dims[1]; Z_DIM = dims[2]

	# apply a gaussian filter
	all_array = ndimage.gaussian_filter(all_array, sigma)

	# threshold to a binary bead mask
	# otsu filter https://en.wikipedia.org/wiki/Otsu%27s_method
	if threshold in ('per-slice', 'per-slice-truncated'):
		num_slice = all_array.shape[2]
		if threshold == 'per-slice-truncated':
			num_slice = signal_depth(all_array)
		# apply an otsu filter, specify the filter at each z slice. Slices at
		# and above num_slice are left zero, which is what zeroing their
		# voxels would produce and avoids running Otsu on an empty slice.
		bw = np.zeros((all_array.shape))
		for kk in range(0,num_slice):
			thresh = threshold_otsu(all_array[:,:,kk])
			bw[:,:,kk] = all_array[:,:,kk] > thresh
	else:
		if threshold == 'global':
			thresh = threshold_otsu(all_array)
		elif isinstance(threshold, str):
			raise ValueError(
				"threshold must be 'per-slice', 'global', or a number; got %r"
				% (threshold,))
		else:
			thresh = float(threshold)
		# Kept as float, matching the per-slice branch's dtype, so `label` and
		# `regionprops` see the same input either way.
		bw = (all_array > thresh).astype(float)

	# find connected volumes within the image, assume each connected volume is a bead
	# record the centroid of each connected volume as the location of the beads
	# relies on https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops
	label_img = label(bw, connectivity=bw.ndim)
	props = regionprops(label_img)
	centroids = np.zeros((len(props),3))
	for kk in range(len(props)):
		centroids[kk]=props[kk].centroid

	centroids_order = np.zeros(centroids.shape)
	centroids_order[:,0] = centroids[:,1] * X_DIM / bw.shape[1]
	centroids_order[:,1] = centroids[:,0] * Y_DIM / bw.shape[0]
	centroids_order[:,2] = centroids[:,2] * Z_DIM / bw.shape[2]

	beads = fmbeads.FMBeads(points=centroids_order)

	return beads

