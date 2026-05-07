import os
import re
from langchain_core.prompts import PromptTemplate

def load_prompt_template(prompt_name: str):
    """
    Loads a prompt template from a Markdown file and processes it for use with LangChain's PromptTemplate.
    
    Args:
        prompt_name: The name of the prompt template file (without extension).    
    Returns:
        system_prompt: The processed template string, with placeholders converted for LangChain use.
        input_vars: The input variables name in the prompt.
    """


    # Construct the path to the template file using the script's directory and prompt name
    template_path = os.path.join(os.path.dirname(__file__), f"{prompt_name}.md")
    
    try:
        # Read the content of the file with proper file handling and UTF-8 encoding
        with open(template_path, 'r', encoding='utf-8') as file:
            template = file.read()
    except FileNotFoundError:
        # Raise an informative error if the file is not found
        raise FileNotFoundError(f"Prompt template '{prompt_name}.md' not found in {os.path.dirname(__file__)}")
    except Exception as e:
        # Catch any other exceptions during file reading, such as permission errors, and raise a ValueError
        raise ValueError(f"Error reading prompt template: {e}")
    

    # # Find all input variables marked with <<variable>> using regex
    input_vars = re.findall(r"<<([^>>]+)>>", template)
    # Replace { / } with {{ / }} for templating compatibility
    template = template.replace("{", "{{").replace("}", "}}")
    # Replace << / >> with { / } for templating input variables
    template = template.replace("<<", "{").replace(">>", "}")
    
    # Convert to system prompt using langchain
    system_prompt = PromptTemplate(
        input_variables = input_vars,
        template = template,
        )

    return system_prompt, input_vars