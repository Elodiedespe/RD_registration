import os
import subprocess

# Affine registration for 5_9 yers old atlas

atlas = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/ANTS9-5Years3T_brain.nii.gz"
atlas_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/ANTS9-5Years3T_head.nii.gz"
atlas_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/atlas_regis_head.nii"
trans_aff = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/atlas_regis_head.txt"
label = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/ANTS9-5Years3T_brain_ANTS_LPBA40_atlas.nii.gz"
label_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_5_9/ANTS9-5Years3T_label_regis_head.nii"


cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", atlas,
        "-ref", atlas_head, "-out", atlas_regis_head]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)

cmd = ["flirt", "-in", label,
        "-ref", atlas_head, "-applyxfm", "-init", trans_aff, "-out", label_regis_head, "-interp", "nearestneighbour"]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)


# Affine registration for 2-5 yers old atlas

atlas = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/ANTS2-5Years_brain.nii.gz"
atlas_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/ANTS2-5Years_head.nii.gz"
atlas_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/ANTS2-5Years_regis_head.nii"
trans_aff = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/atlas_regis_head.txt"
label = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/ANTS2-5Years_brain_LPBA40_atlas.nii.gz"
label_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_2_5/ANTS2-5Years_label_regis_head.nii"

cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", atlas,
        "-ref", atlas_head, "-out", atlas_regis_head]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)

cmd = ["flirt", "-in", label,
        "-ref", atlas_head, "-applyxfm", "-init", trans_aff, "-out", label_regis_head, "-interp", "nearestneighbour"]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)

# Affine registration for 0-2 yers old atlas

atlas = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/ANTS2-0Years_brain.nii.gz"
atlas_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/ANTS2-0Years_head.nii.gz"
atlas_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/atlas_regis_head.nii"
trans_aff = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/atlas_regis_head.txt"
label = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/ANTS2-0Years_brain_ANTS_LPBA40_atlas.nii.gz"
label_regis_head = "/media/mfpgt/MainDataDisk/Elodie/ants_test/test_registration_atlas/atlas/atlas_0-2/ANTS2-0Years_label_regis_head.nii"

cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", atlas,
        "-ref", atlas_head, "-out", atlas_regis_head]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)

cmd = ["flirt", "-in", label,
        "-ref", atlas_head, "-applyxfm", "-init", trans_aff, "-out", label_regis_head, "-interp", "nearestneighbour"]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)
