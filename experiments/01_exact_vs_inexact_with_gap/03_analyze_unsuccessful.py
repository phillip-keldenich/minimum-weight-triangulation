import os, re, bz2
from data_handling import clean_task, classify_runs_in, instance_type_directories


if __name__ == "__main__":
    instance_dirs = instance_type_directories()
    runs_to_delete = []

    delete_classes = {}

    for instance_dir in instance_dirs:
        classified = classify_runs_in(instance_dir)
        for k in classified.keys():
            if k == "UNKNOWN":
                print("UNKNOWN FAILURE MODE FOR RUN:", instance_dir, classified[k][0])
                raise RuntimeError()
            elif k == "GUROBI_UNAVAILABLE" or k == "INVALID_INSTANCE" or k == "TASK_DIED":
                for rname in classified[k]:
                    runs_to_delete.append((instance_dir, rname))
                    if k not in delete_classes: delete_classes[k] = 0
                    delete_classes[k] += 1

    if runs_to_delete:
        print("FOUND", len(runs_to_delete), "RUNS TO DELETE...")
        for delete_class, count in delete_classes.items():
            print("\t", delete_class, ":", count, " runs")
        answer = input("Continue? [y/n]")
        if answer.lower() == 'y' or answer.lower() == 'z':
            for idir, rname in runs_to_delete:
                clean_task(idir, rname)
        else:
            print("Aborting!")
    else:
        print("No runs to delete!")
