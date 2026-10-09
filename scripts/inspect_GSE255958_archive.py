import tarfile

archive_path = "data/raw/GSE255958/GSE255958_RAW.tar"

with tarfile.open(archive_path, "r") as tar:
    members = tar.getmembers()

    print("Number of files:", len(members))
    print("\nFiles inside archive:\n")

    for member in members:
        print(member.name)