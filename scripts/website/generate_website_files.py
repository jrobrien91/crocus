import os
import json

def scan_directory(directory):
    """ Recursively scans a directory and returns a nested dictionary of files. """
    directory_structure = {}

    for entry in os.listdir(directory):
        entry_path = os.path.join(directory, entry)

        if os.path.isdir(entry_path):
            # Recursively scan subdirectories
            directory_structure[entry] = scan_directory(entry_path)
        elif os.path.isfile(entry_path):
            # Add files to the list
            directory_structure.setdefault("", []).append(entry)

    return directory_structure

# Define the root folders to scan
folders = ["CROCUS"]  # Add other root folders if needed
file_structure = {}

for folder in folders:
    if os.path.exists(folder):
        file_structure[folder] = scan_directory(folder)

# Save to files.json
with open("files.json", "w") as f:
    json.dump(file_structure, f, indent=2)

print("✅ files.json updated with subfolders and files!")

