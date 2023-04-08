
for folder in wheelhouse/*; do
    # check if folder is a directory
    if [[ -d "$folder" ]]; then
        # move all files and folders inside folder to parent directory
        echo  moving "$folder"/*  to dist/
        cp -r "$folder"/* dist/
        # remove empty folder
        # rmdir "$folder"
    fi
done