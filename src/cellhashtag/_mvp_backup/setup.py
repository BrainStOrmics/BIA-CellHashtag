import os 
import subprocess
from pathlib import Path

from langchain_openai import ChatOpenAI
from ..config import *

current_file_path = Path(os.path.abspath(__file__))

def setup_LLMs() -> None:
    """
    Initialize and configure all Large Language Models for the BIA-Ghostcoder system.
    
    This function sets up three types of models essential for bioinformatics analysis:
    - Chat model: For natural language understanding and task interpretation
    - Code model: For generating bioinformatics analysis code
    - Embedding model: For vector similarity search in code retrieval
    
    The function handles multiple embedding providers (DashScope, OpenAI) with fallback
    mechanisms to ensure robust model initialization even if some services are unavailable.
    
    Returns:
        None: Models are stored in the global llm_config.MODELS dictionary.
    """
    # Initialize chat model for natural language processing and task understanding
    try:
        chat_model = ChatOpenAI(
            api_key=llm_config.CHAT_MODEL_API['api'],
            base_url=llm_config.CHAT_MODEL_API['url'],
            model=llm_config.CHAT_MODEL_API['model'],
            temperature=0,  # Deterministic output for consistent behavior
            max_retries=3,  # Robust error handling for network issues
            extra_body = {
                "enable_thinking": llm_config.ENABLE_THINKING,
                "enable_search": llm_config.ENABLE_SEARCH,
            }
        )
    except Exception as e:
        # Graceful degradation - system can still function without chat model
        print(f"Warning: Failed to initialize chat model: {e}")
        chat_model = None

    # Store chat model in global configuration for system-wide access
    llm_config.MODELS['chat_model'] = chat_model

    # Initialize multi-modal chat model with base_64 pic imput
    try:
        mmchat_model = ChatOpenAI(
            api_key=llm_config.MULTIMODAL_CHAT_MODEL_API['api'],
            base_url=llm_config.MULTIMODAL_CHAT_MODEL_API['url'],
            model=llm_config.MULTIMODAL_CHAT_MODEL_API['model'],
            temperature=0,  # Deterministic code generation for reproducibility
            max_retries=3,  # Robust error handling for network issues
            extra_body = {
                "enable_thinking": llm_config.ENABLE_THINKING,
                "enable_search": llm_config.ENABLE_SEARCH,
            }
        )
    except Exception as e:
        # Graceful degradation - system can still function without code model
        print(f"Warning: Failed to initialize Multimoda model: {e}")
        mmchat_model = None

    # Store code model in global configuration
    llm_config.MODELS['mmchat_model'] = mmchat_model

def setup_brief(
        config_yaml = 'DEFAULT'# Use default or abs path of your .ymal
        ):
    # Load config
    if config_yaml == 'DEFAULT':
        load_yaml_config(
            os.path.join(current_file_path.parents[1],"config/config.yaml")
        )
        
    else: 
        try:
            load_yaml_config(config_yaml)
        except FileNotFoundError as e:
            print(f"Error loading config file: {e}")
            raise e 

    # Set up LLMs
    setup_LLMs()

    print("\n" + "="*50)
    print("BIA-Brief System Initialization Summary")
    print("="*50)

    # 1. LLM model status
    print("\n[LLM MODELS]")
    print(f"  Chat Model: {'✅ Initialized' if llm_config.MODELS['chat_model'] else '❌ Failed'}")
    print(f"  Multimoda Chat Model: {'✅ Initialized' if llm_config.MODELS['mmchat_model'] else '❌ Failed'}")
    
    print("\n" + "="*50 + "\n")
    print("System initialization complete..")