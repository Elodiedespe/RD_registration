def get_roi_mask(roi_filename, label_number=1):
    mask = nib.load(roi_filename).get_data()
    mask = mask.reshape(mask.shape[0:3])
    mask = mask == label_number
    return mask


def extract_data_in_mask(imgs, mask):
    """ Intersect the imgs with the specified ROI

    The current code assumes that mask and the imgs have the same shape!

    Parameters
    ----------
    imgs : is a list of 3D img files,
    roi_filename : is the name of an img file containing the ROI mask

    Returns
    -------
    data : array with imgs information in the given ROI
    """
    if type(imgs) is not list:
        imgs = [imgs]

    data = np.zeros([mask.shape[0], mask.shape[1], mask.shape[2], len(imgs)])
    for i in range(len(imgs)):
        assert os.path.isfile(imgs[i])
        img = nib.load(imgs[i]).get_data()
        img = img.reshape(img.shape[:3])
        assert img.shape == mask.shape
        data[:, :, :, i] = img

    return data[mask, :]

# custom way
mask = get_roi_mask(path, label)
data = extract_data_in_mask(radio_img, mask)
np.max(data)
np.mean(data)

# nilearn way
masker = NiftiMasker(mask)
data = masker.fit_transform()
result_img = masker.inverse_transform(results)
nib.save(result_img, path)
