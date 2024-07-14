#!/usr/bin/env python3
import os
import sys
import json

def main():
    if len(sys.argv) != 2:
        print("Usage: script.py <path_to_directory>")
        sys.exit(1)
    
    directory = sys.argv[1]
    
    # Walk through all the subdirectories in the given directory
    for root, dirs, files in os.walk(directory):
        if '.cargo-checksum.json' in files:
            json_file_path = os.path.join(root, '.cargo-checksum.json')
            with open(json_file_path, 'r') as file:
                data = json.load(file)

            # Check if 'files' key is in data
            if 'files' in data:
                original_files = data['files']
                filtered_files = {
                    key: value for key, value in original_files.items()
                    if not (key.endswith('.md') or 'tests/' in key or 'test/' in key or 'examples/' in key or 'benches/' in key)
                }
                data['files'] = filtered_files

                # Write the modified data back to the JSON file
                with open(json_file_path, 'w') as file:
                    json.dump(data, file)

if __name__ == "__main__":
    main()

