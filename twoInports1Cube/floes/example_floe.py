from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileInputCube
from cubes.cubes import TestCube 

job = WorkFloe("Example")

job.description = """
**Title Goes Here**

And maybe some description

"""

job.classification = [
    ["OpenEye", "Example"],
    ["OpenEye", "Custom"]
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
tfs = FileInputCube("tfs" )
twoInportsCube = TestCube( "TwoInports" )

# Promotes the parameter
ifs.promote_parameter("data_in", promoted_name="ifs")
tfs.promote_parameter("name", promoted_name="tfs")

job.add_cubes(ifs, tfs, twoInportsCube)
#job.add_cubes(ifs, twoInportsCube)

ifs.success.connect(twoInportsCube.intake)
tfs.success.connect(twoInportsCube.txtintake)

if __name__ == "__main__":
    job.run()
