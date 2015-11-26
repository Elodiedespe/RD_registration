
# System import
import nibabel
import numpy

ct = nibabel.load("ct.nii.gz")
rd = nibabel.load("rdm.nii.gz")
rd_data = rd.get_data()
cta = ct.get_affine()
rda = rd.get_affine()

def inverse_affine(affine):
    """ Invert an affine transformation.
    """
    invr = numpy.linalg.inv(affine[:3, :3])
    inv_affine = numpy.zeros((4, 4))
    inv_affine[3, 3] = 1
    inv_affine[:3, :3] = invr
    inv_affine[:3, 3] =  - numpy.dot(invr, affine[:3, 3])
    return inv_affine

icta = inverse_affine(cta)
irda = inverse_affine(rda)

t = numpy.dot(irda, cta)

print numpy.dot(t, numpy.array([260, 268, 205, 1]))
#print stop

print ct.shape
rd_rescale = numpy.zeros(ct.shape)
for x in range(ct.shape[0]):
    for y in range(ct.shape[1]):
        for z in range(ct.shape[2]):
            voxel_ct = numpy.array([x, y, z, 1])
            voxel_rd = numpy.dot(t, voxel_ct)[:3]
            if (voxel_rd > 0).all() and (voxel_rd < (numpy.asarray(rd.shape) - 1)).all():
                rd_voxel = numpy.round(voxel_rd)
                rd_rescale[x, y, z] = rd_data[rd_voxel[0], rd_voxel[1], rd_voxel[2]]

rd_rescale_im = nibabel.Nifti1Image(rd_rescale, cta)
nibabel.save(rd_rescale_im, "rd_test.nii.gz")





