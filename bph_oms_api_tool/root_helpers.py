import ROOT
def list_directories(root_file, parent_dir="", depth=0):
    # Print the current directory
    print("  " * depth + parent_dir)

    # Get the keys in the current directory
    keys = root_file.GetListOfKeys()

    # Iterate over keys
    for key in keys:
        obj = key.ReadObj()

        # Check if the object is a directory
        if isinstance(obj, ROOT.TDirectory):
            # Recursively call list_directories for subdirectories
            list_directories(obj, parent_dir + "/" + key.GetName(), depth + 1)


 
