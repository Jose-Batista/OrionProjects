from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from chunkdb_modules.cubes import SleepyCube

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
ofs = OEMolOStreamCube("ofs")

# Promotes the parameter
ifs.promote_parameter("data_in", promoted_name="ifs")

# Makes the value static
ofs.set_parameters(data_out="out.csv")

sleeper = SleepyCube("sleep")

job.add_cubes(ifs, ofs, sleeper)
ifs.success.connect(sleeper.intake)
sleeper.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
