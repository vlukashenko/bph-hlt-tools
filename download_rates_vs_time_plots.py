import json 

import requests

def download_image(url, filename):
    """
    Downloads an image from a given URL and saves it to a local file.

    Args:
    url (str): The URL of the image.
    filename (str): The name of the file where the image will be saved.

    Returns:
    bool: True if download was successful, False otherwise.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raises an exception if the request returned an error
        with open(filename, 'wb') as f:
            f.write(response.content)
        return True
    except Exception as e:
        print(f"An error occurred: {e}")
        return False



import re

def remove_version_number(s):
    """
    Removes the integer following '_v' in the provided string.

    Args:
    s (str): The input string containing the '_v' followed by an integer.

    Returns:
    str: A modified string with the integer removed after '_v'.
    """
    # Pattern to find '_v' followed by one or more digits
    pattern = r'(_v)\d+'
    # Replace the matched pattern with '_v'
    result = re.sub(pattern, r'\1', s)
    return result

with open('Path_filters_Menu24_B.json', 'r') as input:
     lines = json.load(input)

names = list(lines.keys())

standard_url = "https://sdonato.web.cern.ch/OMSRatesNtuple/OMSRatesNtuple/RatePlots/plots/2024_physics_allHLT/"

for n in names:
    trigger = remove_version_number(n)
    download_image(standard_url+"rates_"+trigger+"_vsTime.png", "plots/rates_"+trigger+"_vsTime.png")


