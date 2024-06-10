from pathlib import Path
import os
import re
import sys
import io
import json 


def findFile(directory, pattern):
    compiled_pattern = re.compile(pattern)
    for filename in os.listdir(directory):
        if compiled_pattern.match(filename):
            return os.path.join(directory, filename)
    return None     

def expandPath(path):
    return os.path.join(os.getcwd(), path)

def captureObject(obj):
    buffer = io.StringIO()
    original_stdout = sys.stdout
    sys.stdout = buffer
    print(obj)
    sys.stdout = original_stdout
    captured_output = buffer.getvalue()
    return captured_output


def harvestInfoJSON(file_name, pattern_to_harvest):
    if not file_name.endswith('.json'):
        raise ValueError("The file must be in JSON format with a .json extension")
    
    with open(file_name, 'r') as file:
        data = json.load(file)
    
    def find_pattern(data, pattern):
        if isinstance(data, dict):
            for key, value in data.items():
                if key == pattern:
                    yield value
                elif isinstance(value, dict):
                    yield from find_pattern(value, pattern)
                elif isinstance(value, list):
                    for item in value:
                        yield from find_pattern(item, pattern)

    results = set()
    for value in find_pattern(data, pattern_to_harvest):
        parts = value.replace("OR", ";").replace("AND", ";").split(";")
        cleaned_parts = [part.strip() for part in parts if part.strip()]
        results.update(cleaned_parts)

    return list(results)


def flattenList(nested_list):
    stack = nested_list[::-1]
    flattened_list = []
    while stack:
        item = stack.pop()
        if isinstance(item, list):
            stack.extend(item[::-1])
        else:
            flattened_list.append(item)
    return flattened_list


def removeDuplicates(input_list):
    seen = set()
    return [x for x in input_list if not (x in seen or seen.add(x))]
