import platform
import shutil
import subprocess

with open(snakemake.output.unsorted, "wb") as wfd:
    for f in [
        snakemake.input["areashydro"],
        snakemake.input["areassolar"],
        snakemake.input["areaswindon"],
        snakemake.input["areaswindoff"],
        snakemake.input["areaswindofffloating"],
    ]:
        with open(f, "rb") as fd:
            shutil.copyfileobj(fd, wfd)

if platform.system() == "Linux":
    print("linux")
    subprocess.run(
        ["sort", "-g", "-o", snakemake.output.areassorted, snakemake.output.unsorted],
        check=True,
    )

if platform.system() == "Windows":
    print("windows")
    subprocess.run(
        ["sort", "/o", snakemake.output.areassorted, snakemake.output.unsorted],
        check=True,
    )

with open(snakemake.output[0], "wb") as wfd:
    for f in [
        snakemake.input["vreareaheader"],
        snakemake.output.areassorted,
        snakemake.input["genericfooter"],
    ]:
        with open(f, "rb") as fd:
            shutil.copyfileobj(fd, wfd)
