import dataclasses
import enum
import itertools
import os
import shutil
import sys
import typer
from pathlib import Path
from typing import List

from hyperqueue import Client, Job
from hyperqueue.job import ResourceRequest

ROOT = "HyperQueueScriptsCPU"

def run_PS(
        client: Client,
        data_directory: str,
        simulations: int,
):
    job = Job(ROOT + "/../HyperQueueRuns/startHQ/hq-servers")

    shutil.rmtree(data_directory + "/hq-work", ignore_errors=True)
    os.mkdir(data_directory + "/hq-work")
    cwd = data_directory + "/hq-work/job-%{JOB_ID}/"
    job.program(
        args=[ROOT + "/hq-md-compute.sh"],
        cwd=cwd,
        stdout="%{CWD}/stdout",
        stderr="%{CWD}/stderr",
        env=dict(CWD=data_directory,
                 NUM_ITERS=str(simulations)
                 ),
        resources=ResourceRequest(cpus=128)
    )

    # Submit the job
    submitted = client.submit(job)

    # Wait until the job completes
    client.wait_for_jobs([submitted])


def main(data_directory: Path):
    client = Client(ROOT + "/../HyperQueueRuns/startHQ/hq-servers")

    total_sims = 8

    # We might have to restart from partially completed and analysed directories
    # Exclude complete poses and analysis files
    poses = [f for f in os.listdir(data_directory) if not f.endswith("tar.gz") and not f.startswith("failed")]

    for i in range(len(poses)):
        try:
            run_PS(
                client,
                str(data_directory) + "/" + poses[i],
                total_sims
            )
        except:
            print("Interrupting PS workflow run")
            pass


if __name__ == "__main__":
    typer.run(main)
