import base64
from pathlib import Path



def check_file_exists(image_path):
    try:
        with open(image_path, "rb"):
            return True
    except FileNotFoundError:
        return False

# local image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        image_base64 = base64.b64encode(image_file.read()).decode("utf-8")
    return image_base64


# 帮我写一个代码，以文本形式读取代码文件
def read_code_file(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        code_content = file.read()
    return code_content

def tree_dir(
        target_dir, 
        prefix: str = '',
        exclude_dirs: list[str] = None):
    "print directory tree"
    if exclude_dirs is None:
        # Exclude unused directories
        exclude_dirs = ['__pycache__', 'venv', '.venv']
        
    path = Path(target_dir)
    contents = sorted([p for p in path.iterdir() if p.name not in exclude_dirs])
    
    pointers = ['├── '] * (len(contents) - 1) + ['└── ']
    for pointer, child in zip(pointers, contents):
        print(prefix + pointer + child.name)
        if child.is_dir():
            extension = '│   ' if pointer == '├── ' else '    '
            tree_dir(child, prefix + extension, exclude_dirs)

###################
# Old Functions
###################

def parse_jsonfromcontent(content: str):
    '''
    Extracts and parses JSON data from a string containing JSON data.
    This function looks for '```json'' and '```'' tags in the input string and extracts the JSON data between them.
    It then tries to parse the extracted string into Python objects (usually a dictionary or list).

    Parameters: 
    content (str): A string containing the JSON data.

    Returns: 
    dict or list: the parsed Python object.

    Raises: 
    ValueError: If the '```json'' or '```' token is not found, or JSON decoding fails
    '''
    start_marker = '```json'
    end_marker = '```'
    # search JSON start position
    start = content.find(start_marker)
    if start == -1:
        raise ValueError("Cannot find JSON start string: '```json'")
    # calculate JSON data start position
    start += len(start_marker)
    # search JSON end position
    end = content.find(end_marker, start)
    if end == -1:
        raise ValueError("Cannot find JSON end string '```'")
    # extract JSON string and remove whitespace     
    json_str = content[start:end].strip()
    # trying to parse JSON data
    try:
        return json.loads(json_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Failed to parse JSON: {e}") 


