# **Welcome to Cell#**  
Anyone can annotation like a pro with Cell#!


## What is Cell#

Cell# is an automated annotation tool based on Large Language Models (LLMs), specifically designed for single-cell RNA sequencing (scRNAseq) and spatial RNA sequencing (space-RNAseq) data. With the use of self-critic mechanisms and web scraping technologies, Cell# achieves expert-level annotation quality.

### Features

- **Advanced Annotation Capabilities**: Utilizes self-critic and web scraping techniques to achieve high-precision annotations for bioinformatics data.
- **Powerful Framework Support**: Built primarily on the [Langgraph](https://github.com/langgraph) framework, making full use of its flow control and map-reduce functionalities.
- **Flexible Data Compatibility**: Currently supports input data in scanpy and AnnData formats, aiming to help users save time and improve work efficiency.

### Agent Framework
Cell# is an open-source Agent based on community framework [LangGraph](https://langchain-ai.github.io/langgraph/) and LangChain. We applied self-critic (refection), map-reduce and more techs in Cell#. 
![Cognitive architecture of Cell# agent](https://lh3.googleusercontent.com/pw/AP1GczO2yNqxJgqSsSuKpPDjHfyV86wmOab3WRrjmgpkfARFZE5U0EI1kG7803CqC1kETEca4AJiS6vgIKiMOl2eP7n3kN_RX5EjBRNYvh7fHQ-lWlIE1f1GTL29SgkEQ3G-Rg2bi0fHdml81Q8B3SGqZ_5j=w1167-h1115-s-no-gm?authuser=0)

[Join Langchain community here](https://www.langchain.com/join-community) (Hi! I am not invited yet!)
![LangGraph](https://langchain-ai.github.io/langgraph/static/wordmark_dark.svg)

PS: Cell# is an Agent exercise, a component of the BIAgent upgrade program, mainly to provide demos for partners, and to provide code wrapping standards internally. But, I'll be happy to listen to your input to maintain the tool for as long as possible. ♥(´∀` 

## Play with Cell#
 ### Prerequisites
1. Python 3.11+ (I only test on 3.11 # (˘•ω•˘))
2. A Tavily API key (required for web scraping functionality)
3. And I am not yet upload it to pip so... you have to clone this repo to your local

#### Step 1: Clone the Repository
```bash
git clone https://github.com/your-username/cellsharp.git
cd cellsharp
```
#### Step 2: Create a Virtual Environment (Optional but Recommended)
```
# Create a conda environment (adjust Python version as needed)
conda create --name cellhastag_env python=3.11

# Activate the environment
conda activate cellhastag_env 
```
#### Step 3: Install Dependencies
```
# Install core dependencies
pip install scanpy pandas langchain langgraph
```

#### Step 4: Configure Tavily API Key (Optional but Recommended)

1.  Visit  [Tavily's website](https://tavily.com/)  and sign up for free credits.
2.  Save your API key in a file named  `.env`  in the project root:
  ```
  TAVILY_API_KEY=your_api_key_here
  ```

### Cell# usage
Please have a look at the **DEMO.ipynb** Notebook, super easy and... easy

Note:
For people from BGI: 杭州的地区的可以找刘石平老师要api, 其他老师同学可能要去找智能所的老师们要啦。


## Todo list

 - [ ]  **Add Pip Installation Support**
    - [ ] Package Cell# for PyPI to enable installation via  `pip install cellsharp`.
    - [ ] Create documentation for seamless distribution and version management.
- [ ] **Integrate Human-in-the-Loop (HITL) Component**
	- [ ] Allow users to review and refine annotations generated by the LLM.
	- [ ] Implement a feedback loop to improve model performance over time based on user input.
 - [ ] **Expand Web Scraper Support**
	 - [ ] Add support for additional web scraping sources beyond Tavily (e.g., other APIs or databases).
	 - [ ] Improve reliability and error handling for web scraping workflows.
 - [ ] **Normalize Cell Marker Input**
	- [ ] Standardize input formats for cell markers to ensure consistency and reduce user configuration overhead.
	- [ ] Implement validation checks for cell marker data to prevent errors during annotation.

### acknowledgements
![BGI_research](https://research.genomics.cn/public/img/logo-blue.png?v=1) ![enter image description here](https://upload.wikimedia.org/wikipedia/en/thumb/6/64/University_of_the_Chinese_Academy_of_Sciences_logo.svg/220px-University_of_the_Chinese_Academy_of_Sciences_logo.svg.png)


> Written with [StackEdit](https://stackedit.io/).
