import os
import argparse

parser = argparse.ArgumentParser(
                    prog='PolyA_comparison',
                    description='Execute circRNA-sponging pipeline in parallel for totalRNA and reference PolyA data to determine false positive circRNA detection rates')
parser.add_argument('-p', '--polyA', help='circRNA-sponging config file for polyA input data')
parser.add_argument('-t', '--totalRNA', help='circRNA-sponging config file for totalRNA input data')
parser.add_argument('-w', '--workdir', help='Nextflow work directory path')
parser.add_argument('-profile', '--profile', help='Nextflow -profile options, e.g. singularity,cluster')
parser.add_argument('-d', '--dummy', help='Only echo the commands', action='store_true')
args = parser.parse_args()


def call_pipeline(configs_path, profile_opts, workdir, dummy=False):
    cmd = ["nextflow", "run", "./", "-c", configs_path, "-profile", profile_opts, "-w", workdir, "--circRNA_only", "true"]
    print("[info] calling process:", " ".join(cmd))
    os.spawnl(os.P_DETACH, " ".join(cmd))


# read CLI and start pipelines
call_pipeline(args.polyA, args.profile, args.workdir, dummy=args.dummy)
call_pipeline(args.totalRNA, args.profile, args.workdir, dummy=args.dummy)
