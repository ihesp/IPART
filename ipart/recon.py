
import numpy as np


def rank_order(image):
    """Return an image of the same shape where each pixel is the
    index of the pixel value in the ascending order of the unique
    values of `image`, aka the rank-order value.

    Parameters
    ----------
    image: ndarray

    Returns
    -------
    labels: ndarray of type np.uint32, of shape image.shape
        New array where each pixel has the rank-order value of the
        corresponding pixel in `image`. Pixel values are between 0 and
        n - 1, where n is the number of distinct unique values in
        `image`.

    original_values: 1-D ndarray
        Unique original values of `image`

    Examples
    --------
    >>> a = np.array([[1, 4, 5], [4, 4, 1], [5, 1, 1]])
    >>> a
    array([[1, 4, 5],
           [4, 4, 1],
           [5, 1, 1]])
    >>> rank_order(a)
    (array([[0, 1, 2],
           [1, 1, 0],
           [2, 0, 0]], dtype=uint32), array([1, 4, 5]))
    >>> b = np.array([-1., 2.5, 3.1, 2.5])
    >>> rank_order(b)
    (array([0, 1, 2, 1], dtype=uint32), array([-1. ,  2.5,  3.1]))
    """
    flat_image = image.ravel()
    sort_order = flat_image.argsort().astype(np.uint32)
    flat_image = flat_image[sort_order]
    sort_rank = np.zeros_like(sort_order)
    is_different = flat_image[:-1] != flat_image[1:]
    np.cumsum(is_different, out=sort_rank[1:])
    original_values = np.zeros((sort_rank[-1] + 1,), image.dtype)
    original_values[0] = flat_image[0]
    original_values[1:] = flat_image[1:][is_different]
    int_image = np.zeros_like(sort_order)
    int_image[sort_order] = sort_rank
    return (int_image.reshape(image.shape), original_values)

def reconstructionLSM(seed, mask, lsm, landocean, method='dilation', selem=None,
        offset=None):
    """Perform a morphological reconstruction of an image.

    Morphological reconstruction by dilation is similar to basic morphological
    dilation: high-intensity values will replace nearby low-intensity values.
    The basic dilation operator, however, uses a structuring element to
    determine how far a value in the input image can spread. In contrast,
    reconstruction uses two images: a "seed" image, which specifies the values
    that spread, and a "mask" image, which gives the maximum allowed value at
    each pixel. The mask image, like the structuring element, limits the spread
    of high-intensity values. Reconstruction by erosion is simply the inverse:
    low-intensity values spread from the seed image and are limited by the mask
    image, which represents the minimum allowed value.

    Alternatively, you can think of reconstruction as a way to isolate the
    connected regions of an image. For dilation, reconstruction connects
    regions marked by local maxima in the seed image: neighboring pixels
    less-than-or-equal-to those seeds are connected to the seeded region.
    Local maxima with values larger than the seed image will get truncated to
    the seed value.

    Parameters
    ----------
    seed : ndarray
        The seed image (a.k.a. marker image), which specifies the values that
        are dilated or eroded.
    mask : ndarray
        The maximum (dilation) / minimum (erosion) allowed value at each pixel.
    method : {'dilation'|'erosion'}
        Perform reconstruction by dilation or erosion. In dilation (or
        erosion), the seed image is dilated (or eroded) until limited by the
        mask image. For dilation, each seed value must be less than or equal
        to the corresponding mask value; for erosion, the reverse is true.
    selem : ndarray
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    reconstructed : ndarray
       The result of morphological reconstruction.
    """

    assert tuple(seed.shape) == tuple(mask.shape)
    assert tuple(seed.shape) == tuple(lsm.shape)
    if method == 'dilation' and np.any(seed > mask):
        raise ValueError("Intensity of seed image must be less than that "
                         "of the mask image for reconstruction by dilation.")
    elif method == 'erosion' and np.any(seed < mask):
        raise ValueError("Intensity of seed image must be greater than that "
                         "of the mask image for reconstruction by erosion.")
    try:
        from _greyreconstruct import reconstruction_loop
    except ImportError:
        raise ImportError("_greyreconstruct extension not available.")

    if selem is None:
        selem = np.ones([3] * seed.ndim, dtype=bool)
    else:
        selem = selem.astype(bool)

    if offset is None:
        if not all([d % 2 == 1 for d in selem.shape]):
            raise ValueError("Footprint dimensions must all be odd")
        offset = np.array([d // 2 for d in selem.shape])
    # Cross out the center of the selem
    selem[tuple(slice(d, d + 1) for d in offset)] = False

    # Make padding for edges of reconstructed image so we can ignore boundaries
    padding = (np.array(selem.shape) / 2).astype(int)
    dims = np.zeros(seed.ndim + 1, dtype=int)
    dims[1:] = np.array(seed.shape) + 2 * padding
    dims[0] = 2
    inside_slices = tuple(slice(p, -p) for p in padding)
    # Set padded region to minimum image intensity and mask along first axis so
    # we can interleave image and mask pixels when sorting.
    if method == 'dilation':
        pad_value = np.min(seed)
    elif method == 'erosion':
        pad_value = np.max(seed)
    else:
        raise ValueError("Reconstruction method can be one of 'erosion' "
                         "or 'dilation'. Got '%s'." % method)
    images = np.ones(dims) * pad_value
    images[(0, ) + inside_slices] = seed
    images[(1, ) + inside_slices] = mask

    lsmpad = np.ones(dims[1:]) * 0
    lsmpad[inside_slices] = lsm

    # Create a list of strides across the array to get the neighbors within
    # a flattened array
    value_stride = np.array(images.strides[1:]) // images.dtype.itemsize
    image_stride = images.strides[0] // images.dtype.itemsize
    selem_mgrid = np.mgrid[[slice(-o, d - o)
                            for d, o in zip(selem.shape, offset)]]
    selem_offsets = selem_mgrid[:, selem].transpose()
    nb_strides = np.array([np.sum(value_stride * selem_offset)
                           for selem_offset in selem_offsets], np.int32)

    images = images.flatten()
    lsmpad = lsmpad.flatten().astype(np.int32)

    # Erosion goes smallest to largest; dilation goes largest to smallest.
    index_sorted = np.argsort(images).astype(np.int32)
    if method == 'dilation':
        index_sorted = index_sorted[::-1]

    # Make a linked list of pixels sorted by value. -1 is the list terminator.
    prev = -np.ones(len(images), np.int32)
    next = -np.ones(len(images), np.int32)
    prev[index_sorted[1:]] = index_sorted[:-1]
    next[index_sorted[:-1]] = index_sorted[1:]

    # Cython inner-loop compares the rank of pixel values.
    if method == 'dilation':
        value_rank, value_map = rank_order(images)
    elif method == 'erosion':
        value_rank, value_map = rank_order(-images)
        value_map = -value_map

    start = index_sorted[0]
    reconstruction_loop(value_rank, lsmpad, prev, next, nb_strides, start,
                        image_stride, landocean)

    # Reshape reconstructed image to original image shape and remove padding.
    rec_img = value_map[value_rank[:image_stride]]
    rec_img.shape = np.array(seed.shape) + 2 * padding
    return rec_img[inside_slices]

def reconstruction(seed, mask, method='dilation', selem=None, offset=None):
    """Perform a morphological reconstruction of an image.

    Morphological reconstruction by dilation is similar to basic morphological
    dilation: high-intensity values will replace nearby low-intensity values.
    The basic dilation operator, however, uses a structuring element to
    determine how far a value in the input image can spread. In contrast,
    reconstruction uses two images: a "seed" image, which specifies the values
    that spread, and a "mask" image, which gives the maximum allowed value at
    each pixel. The mask image, like the structuring element, limits the spread
    of high-intensity values. Reconstruction by erosion is simply the inverse:
    low-intensity values spread from the seed image and are limited by the mask
    image, which represents the minimum allowed value.

    Alternatively, you can think of reconstruction as a way to isolate the
    connected regions of an image. For dilation, reconstruction connects
    regions marked by local maxima in the seed image: neighboring pixels
    less-than-or-equal-to those seeds are connected to the seeded region.
    Local maxima with values larger than the seed image will get truncated to
    the seed value.

    Parameters
    ----------
    seed : ndarray
        The seed image (a.k.a. marker image), which specifies the values that
        are dilated or eroded.
    mask : ndarray
        The maximum (dilation) / minimum (erosion) allowed value at each pixel.
    method : {'dilation'|'erosion'}
        Perform reconstruction by dilation or erosion. In dilation (or
        erosion), the seed image is dilated (or eroded) until limited by the
        mask image. For dilation, each seed value must be less than or equal
        to the corresponding mask value; for erosion, the reverse is true.
    selem : ndarray
        The neighborhood expressed as an n-D array of 1's and 0's.
        Default is the ball of radius 1 according to the maximum norm
        (i.e. a 3x3 square for 2D images, a 3x3x3 cube for 3D images, etc.)

    Returns
    -------
    reconstructed : ndarray
       The result of morphological reconstruction.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.morphology import reconstruction

    First, we create a sinusoidal mask image with peaks at middle and ends.

    >>> x = np.linspace(0, 4 * np.pi)
    >>> y_mask = np.cos(x)

    Then, we create a seed image initialized to the minimum mask value (for
    reconstruction by dilation, min-intensity values don't spread) and add
    "seeds" to the left and right peak, but at a fraction of peak value (1).

    >>> y_seed = y_mask.min() * np.ones_like(x)
    >>> y_seed[0] = 0.5
    >>> y_seed[-1] = 0
    >>> y_rec = reconstruction(y_seed, y_mask)

    The reconstructed image (or curve, in this case) is exactly the same as the
    mask image, except that the peaks are truncated to 0.5 and 0. The middle
    peak disappears completely: Since there were no seed values in this peak
    region, its reconstructed value is truncated to the surrounding value (-1).

    As a more practical example, we try to extract the bright features of an
    image by subtracting a background image created by reconstruction.

    >>> y, x = np.mgrid[:20:0.5, :20:0.5]
    >>> bumps = np.sin(x) + np.sin(y)

    To create the background image, set the mask image to the original image,
    and the seed image to the original image with an intensity offset, `h`.

    >>> h = 0.3
    >>> seed = bumps - h
    >>> background = reconstruction(seed, bumps)

    The resulting reconstructed image looks exactly like the original image,
    but with the peaks of the bumps cut off. Subtracting this reconstructed
    image from the original image leaves just the peaks of the bumps

    >>> hdome = bumps - background

    This operation is known as the h-dome of the image and leaves features
    of height `h` in the subtracted image.

    Notes
    -----
    The algorithm is taken from [1]_. Applications for greyscale reconstruction
    are discussed in [2]_ and [3]_.

    References
    ----------
    .. [1] Robinson, "Efficient morphological reconstruction: a downhill
           filter", Pattern Recognition Letters 25 (2004) 1759-1767.
    .. [2] Vincent, L., "Morphological Grayscale Reconstruction in Image
           Analysis: Applications and Efficient Algorithms", IEEE Transactions
           on Image Processing (1993)
    .. [3] Soille, P., "Morphological Image Analysis: Principles and
           Applications", Chapter 6, 2nd edition (2003), ISBN 3540429883.
    """
    assert tuple(seed.shape) == tuple(mask.shape)
    if method == 'dilation' and np.any(seed > mask):
        raise ValueError("Intensity of seed image must be less than that "
                         "of the mask image for reconstruction by dilation.")
    elif method == 'erosion' and np.any(seed < mask):
        raise ValueError("Intensity of seed image must be greater than that "
                         "of the mask image for reconstruction by erosion.")
    try:
        #from ._greyreconstruct import reconstruction_loop
        from _greyreconstruct import reconstruction_loop
    except ImportError:
        raise ImportError("_greyreconstruct extension not available.")

    if selem is None:
        selem = np.ones([3] * seed.ndim, dtype=bool)
    else:
        selem = selem.astype(bool)

    if offset is None:
        if not all([d % 2 == 1 for d in selem.shape]):
            raise ValueError("Footprint dimensions must all be odd")
        offset = np.array([d // 2 for d in selem.shape])
    # Cross out the center of the selem
    selem[tuple(slice(d, d + 1) for d in offset)] = False

    # Make padding for edges of reconstructed image so we can ignore boundaries
    padding = (np.array(selem.shape) / 2).astype(int)
    dims = np.zeros(seed.ndim + 1, dtype=int)
    dims[1:] = np.array(seed.shape) + 2 * padding
    dims[0] = 2
    inside_slices = tuple(slice(p, -p) for p in padding)
    # Set padded region to minimum image intensity and mask along first axis so
    # we can interleave image and mask pixels when sorting.
    if method == 'dilation':
        pad_value = np.min(seed)
    elif method == 'erosion':
        pad_value = np.max(seed)
    else:
        raise ValueError("Reconstruction method can be one of 'erosion' "
                         "or 'dilation'. Got '%s'." % method)
    images = np.ones(dims) * pad_value
    images[(0, ) + inside_slices] = seed
    images[(1, ) + inside_slices] = mask
    xpad = tuple(slice(p, -p) for p in padding[:-1])
    images[(0,) + xpad + (slice(0,1),)] = np.array(seed)[...,-1,None]
    images[(0,) + xpad + (slice(-1,None),)] = np.array(seed)[...,0,None]
    images[(1,) + xpad + (slice(0,1),)] = np.array(mask)[...,-1,None]
    images[(1,) + xpad + (slice(-1,None),)] = np.array(mask)[...,0,None]

    #images[0, 1:-1, 0]=seed[:, -1]
    #images[0, 1:-1, -1]=seed[:, 0]
    #images[1, 1:-1, 0]=mask[:, -1]
    #images[1, 1:-1, -1]=mask[:, 0]

    # Create a list of strides across the array to get the neighbors within
    # a flattened array
    value_stride = np.array(images.strides[1:]) // images.dtype.itemsize
    image_stride = images.strides[0] // images.dtype.itemsize
    selem_mgrid = np.mgrid[[slice(-o, d - o)
                            for d, o in zip(selem.shape, offset)]]
    selem_offsets = selem_mgrid[:, selem].transpose()
    nb_strides = np.array([np.sum(value_stride * selem_offset)
                           for selem_offset in selem_offsets], np.int32)

    images = images.flatten()

    # Erosion goes smallest to largest; dilation goes largest to smallest.
    index_sorted = np.argsort(images).astype(np.int32)
    if method == 'dilation':
        index_sorted = index_sorted[::-1]

    # Make a linked list of pixels sorted by value. -1 is the list terminator.
    prev = -np.ones(len(images), np.int32)
    next = -np.ones(len(images), np.int32)
    prev[index_sorted[1:]] = index_sorted[:-1]
    next[index_sorted[:-1]] = index_sorted[1:]

    # Cython inner-loop compares the rank of pixel values.
    if method == 'dilation':
        value_rank, value_map = rank_order(images)
    elif method == 'erosion':
        value_rank, value_map = rank_order(-images)
        value_map = -value_map

    start = index_sorted[0]
    reconstruction_loop(value_rank, prev, next, nb_strides, start,
                        image_stride)

    # Reshape reconstructed image to original image shape and remove padding.
    rec_img = value_map[value_rank[:image_stride]]
    rec_img.shape = np.array(seed.shape) + 2 * padding
    return rec_img[inside_slices]
