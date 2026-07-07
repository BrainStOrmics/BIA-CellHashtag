# Build a deep research agent

[Tutorials](/oss/python/deepagents/data-analysis)[Deep Agents](/oss/python/deepagents/data-analysis)

# Build a deep research agent
Copy page

Build a multi-step web research agent with subagent delegationCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

This guide demonstrates how to build a multi-step web research agent from scratch using [Deep Agents](/oss/python/deepagents). The agent decomposes research questions into focused tasks, delegates them to specialized sub-agents, and synthesizes findings into a comprehensive report.
The agent you build will:

{' ' * (self.list_depth - 1)}- Plan research using a todo list

{' ' * (self.list_depth - 1)}- Delegate focused research tasks to sub-agents with isolated context

{' ' * (self.list_depth - 1)}- Assess search results and plan next steps as you gather information

{' ' * (self.list_depth - 1)}- Synthesize findings with proper citations into a final report

The spawned sub-agents will conduct web searches with Tavily, fetching full webpage content for analysis.

### [​](#key-concepts)Key concepts

This tutorial covers:

{' ' * (self.list_depth - 1)}- [Subagents](/oss/python/deepagents/subagents) for parallel, context-isolated research

{' ' * (self.list_depth - 1)}- Custom [tools](/oss/python/langchain/tools) for web search

{' ' * (self.list_depth - 1)}- Multi-step planning with the [built-in planning tool](/oss/python/deepagents/harness#task-planning)

## [​](#prerequisites)Prerequisites

API keys for:

{' ' * (self.list_depth - 1)}- Anthropic (Claude) or Google (Gemini)

{' ' * (self.list_depth - 1)}- [Tavily](https://www.tavily.com/) for web search (optional - free tier sufficient)

{' ' * (self.list_depth - 1)}- [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-deepagents-deep-research) for tracing (optional)

## [​](#setup)Setup

1[](#)

Create project directory
```
mkdir deep-research-agent
cd deep-research-agent

```
2[](#)

Install dependencies
{' ' * (self.list_depth - 1)}- Claude
{' ' * (self.list_depth - 1)}- Geminipipuv
```
pip install deepagents tavily-python httpx markdownify langchain-anthropic langchain-core

```
pipuv
```
pip install deepagents tavily-python httpx markdownify langchain-google-genai langchain-core

```
3[](#)

Set API keys
{' ' * (self.list_depth - 1)}- Claude
{' ' * (self.list_depth - 1)}- Gemini
```
export ANTHROPIC_API_KEY="your_anthropic_api_key"
export TAVILY_API_KEY="your_tavily_api_key"
export LANGSMITH_API_KEY="your_langsmith_api_key" # Optional

```

```
export GOOGLE_API_KEY="your_google_api_key"
export TAVILY_API_KEY="your_tavily_api_key"
export LANGSMITH_API_KEY="your_langsmith_api_key" # Optional

```

## [​](#build-the-agent)Build the agent

Create `agent.py` in your project directory:
1[](#)

Add toolsAdd the custom search tool. The `tavily_search` tool uses Tavily for URL discovery, then fetches full webpage content so the agent can analyze complete sources instead of summaries.
```
import os
from typing import Annotated, Literal

import httpx
from langchain.tools import InjectedToolArg, tool
from markdownify import markdownify
from tavily import TavilyClient

tavily_client = TavilyClient(api_key=os.environ["TAVILY_API_KEY"])

def fetch_webpage_content(url: str, timeout: float = 10.0) -> str:
 """Fetch webpage and convert HTML to markdown."""
 headers = {
 "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
 }
 try:
 response = httpx.get(url, headers=headers, timeout=timeout)
 response.raise_for_status()
 return markdownify(response.text)
 except Exception as e:
 return f"Error fetching {url}: {e!s}"

@tool(parse_docstring=True)
def tavily_search(
 query: str,
 max_results: Annotated[int, InjectedToolArg] = 1,
 topic: Annotated[
 Literal["general", "news", "finance"], InjectedToolArg
 ] = "general",
) -> str:
 """Search the web for information on a given query.

 Uses Tavily to discover relevant URLs, then fetches and returns full webpage content as markdown.

 Args:
 query: Search query to execute
 max_results: Maximum number of results to return (default: 1)
 topic: Topic filter - 'general', 'news', or 'finance' (default: 'general')

 Returns:
 Formatted search results with full webpage content
 """
 search_results = tavily_client.search(
 query,
 max_results=max_results,
 topic=topic,
 )
 result_texts = []
 for result in search_results.get("results", []):
 url = result["url"]
 title = result["title"]
 content = fetch_webpage_content(url)
 result_texts.append(f"## {title}\n**URL:** {url}\n\n{content}\n---")

 return f"Found {len(result_texts)} result(s) for '{query}':\n\n" + "\n".join(
 result_texts
 )

```
2[](#)

Add promptsAdd the orchestrator workflow and sub-agent prompt templates to `agent.py`:
```
RESEARCH_WORKFLOW_INSTRUCTIONS = """# Research Workflow

Follow this workflow for all research requests:

1. **Plan**: Create a todo list with write_todos to break down the research into focused tasks
2. **Save the request**: Use write_file() to save the user's research question to `/research_request.md`
3. **Research**: Delegate research tasks to sub-agents using the task() tool - ALWAYS use sub-agents for research, never conduct research yourself
4. **Synthesize**: Review all sub-agent findings and consolidate citations (each unique URL gets one number across all findings)
5. **Write Report**: Write a comprehensive final report to `/final_report.md` (see Report Writing Guidelines below)
6. **Verify**: Read `/research_request.md` and confirm you've addressed all aspects with proper citations and structure

## Research Planning Guidelines
- Batch similar research tasks into a single TODO to minimize overhead
- For simple fact-finding questions, use 1 sub-agent
- For comparisons or multi-faceted topics, delegate to multiple parallel sub-agents
- Each sub-agent should research one specific aspect and return findings

## Report Writing Guidelines

When writing the final report to `/final_report.md`, follow these structure patterns:

**For comparisons:**
1. Introduction
2. Overview of topic A
3. Overview of topic B
4. Detailed comparison
5. Conclusion

**For lists/rankings:**
Simply list items with details - no introduction needed:
1. Item 1 with explanation
2. Item 2 with explanation
3. Item 3 with explanation

**For summaries/overviews:**
1. Overview of topic
2. Key concept 1
3. Key concept 2
4. Key concept 3
5. Conclusion

**General guidelines:**
- Use clear section headings (## for sections, ### for subsections)
- Write in paragraph form by default - be text-heavy, not just bullet points
- Do NOT use self-referential language ("I found...", "I researched...")
- Write as a professional report without meta-commentary
- Each section should be comprehensive and detailed
- Use bullet points only when listing is more appropriate than prose

**Citation format:**
- Cite sources inline using [1], [2], [3] format
- Assign each unique URL a single citation number across ALL sub-agent findings
- End report with ### Sources section listing each numbered source
- Number sources sequentially without gaps (1,2,3,4...)
- Format: [1] Source Title: URL (each on separate line for proper list rendering)
- Example:

 Some important finding [1]. Another key insight [2].

 ### Sources
 [1] AI Research Paper: https://example.com/paper
 [2] Industry Analysis: https://example.com/analysis
"""

```
See all 63 lines
```
RESEARCHER_INSTRUCTIONS = """You are a research assistant conducting research on the user's input topic. For context, today's date is {date}.

Your job is to use tools to gather information about the user's input topic.
You can use the tavily_search tool to find resources that can help answer the research question.
You can call it in series or in parallel, your research is conducted in a tool-calling loop.

You have access to the tavily_search tool for conducting web searches.

Think like a human researcher with limited time. Follow these steps:

1. **Read the question carefully** - What specific information does the user need?
2. **Start with broader searches** - Use broad, comprehensive queries first
3. **After each search, pause and assess** - Do I have enough to answer? What's still missing?
4. **Execute narrower searches as you gather information** - Fill in the gaps
5. **Stop when you can answer confidently** - Don't keep searching for perfection

**Tool Call Budgets** (Prevent excessive searching):
- **Simple queries**: Use 2-3 search tool calls maximum
- **Complex queries**: Use up to 5 search tool calls maximum
- **Always stop**: After 5 search tool calls if you cannot find the right sources

**Stop Immediately When**:
- You can answer the user's question comprehensively
- You have 3+ relevant examples/sources for the question
- Your last 2 searches returned similar information

After each search, assess results before continuing: What key information did I find? What's missing? Do I have enough to answer? Should I search more or provide my answer?

When providing your findings back to the orchestrator:

1. **Structure your response**: Organize findings with clear headings and detailed explanations
2. **Cite sources inline**: Use [1], [2], [3] format when referencing information from your searches
3. **Include Sources section**: End with ### Sources listing each numbered source with title and URL

Example:
## Key Findings

Context engineering is a critical technique for AI agents [1]. Studies show that proper context management can improve performance by 40% [2].

### Sources
[1] Context Engineering Guide: https://example.com/context-guide
[2] AI Performance Study: https://example.com/study

The orchestrator will consolidate citations from all sub-agents into the final report.
"""

```
See all 45 lines
```
SUBAGENT_DELEGATION_INSTRUCTIONS = """# Sub-Agent Research Coordination

Your role is to coordinate research by delegating tasks from your TODO list to specialized research sub-agents.

## Delegation Strategy

**DEFAULT: Start with 1 sub-agent** for most queries:
- "What is quantum computing?" -> 1 sub-agent (general overview)
- "List the top 10 coffee shops in San Francisco" -> 1 sub-agent
- "Summarize the history of the internet" -> 1 sub-agent
- "Research context engineering for AI agents" -> 1 sub-agent (covers all aspects)

**ONLY parallelize when the query EXPLICITLY requires comparison or has clearly independent aspects:**

**Explicit comparisons** -> 1 sub-agent per element:
- "Compare OpenAI vs Anthropic vs DeepMind AI safety approaches" -> 3 parallel sub-agents
- "Compare Python vs JavaScript for web development" -> 2 parallel sub-agents

**Clearly separated aspects** -> 1 sub-agent per aspect (use sparingly):
- "Research renewable energy adoption in Europe, Asia, and North America" -> 3 parallel sub-agents (geographic separation)
- Only use this pattern when aspects cannot be covered efficiently by a single comprehensive search

## Key Principles
- **Bias towards single sub-agent**: One comprehensive research task is more token-efficient than multiple narrow ones
- **Avoid premature decomposition**: Don't break "research X" into "research X overview", "research X techniques", "research X applications" - just use 1 sub-agent for all of X
- **Parallelize only for clear comparisons**: Use multiple sub-agents when comparing distinct entities or geographically separated data

## Parallel Execution Limits
- Use at most {max_concurrent_research_units} parallel sub-agents per iteration
- Make multiple task() calls in a single response to enable parallel execution
- Each sub-agent returns findings independently

## Research Limits
- Stop after {max_researcher_iterations} delegation rounds if you haven't found adequate sources
- Stop when you have sufficient information to answer comprehensively
- Bias towards focused research over exhaustive exploration"""

```
See all 36 lines3[](#)

Create the agentAdd the model initialization and agent creation to `agent.py`. Choose your provider:
{' ' * (self.list_depth - 1)}- Claude
{' ' * (self.list_depth - 1)}- GeminiGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from datetime import datetime

from deepagents import create_deep_agent
from langchain.chat_models import init_chat_model

max_concurrent_research_units = 3
max_researcher_iterations = 3

current_date = datetime.now().strftime("%Y-%m-%d")

INSTRUCTIONS = (
 RESEARCH_WORKFLOW_INSTRUCTIONS
 + "\n\n"
 + "=" * 80
 + "\n\n"
 + SUBAGENT_DELEGATION_INSTRUCTIONS.format(
 max_concurrent_research_units=max_concurrent_research_units,
 max_researcher_iterations=max_researcher_iterations,
 )
)

research_sub_agent = {
 "name": "research-agent",
 "description": "Delegate research to the sub-agent. Give one topic at a time.",
 "system_prompt": RESEARCHER_INSTRUCTIONS.format(date=current_date),
 "tools": [tavily_search],
}

model = init_chat_model(model="google_genai:gemini-3.5-flash", temperature=0.0)

agent = create_deep_agent(
 model=model,
 tools=[tavily_search],
 system_prompt=INSTRUCTIONS,
 subagents=[research_sub_agent],
)

```

```
from datetime import datetime

from langchain_google_genai import ChatGoogleGenerativeAI
from deepagents import create_deep_agent

max_concurrent_research_units = 3
max_researcher_iterations = 3

current_date = datetime.now().strftime("%Y-%m-%d")

INSTRUCTIONS = (
 RESEARCH_WORKFLOW_INSTRUCTIONS
 + "\n\n"
 + "=" * 80
 + "\n\n"
 + SUBAGENT_DELEGATION_INSTRUCTIONS.format(
 max_concurrent_research_units=max_concurrent_research_units,
 max_researcher_iterations=max_researcher_iterations,
 )
)

research_sub_agent = {
 "name": "research-agent",
 "description": "Delegate research to the sub-agent. Give one topic at a time.",
 "system_prompt": RESEARCHER_INSTRUCTIONS.format(date=current_date),
 "tools": [tavily_search],
}

model = ChatGoogleGenerativeAI(model="gemini-3-pro-preview", temperature=0.0)

agent = create_deep_agent(
 model=model,
 tools=[tavily_search],
 system_prompt=INSTRUCTIONS,
 subagents=[research_sub_agent],
)

```

## [​](#run-the-agent)Run the agent

You can run the agent synchronously, meaning it will wait for the full result and then print it, or you can stream updates as they come in.
Add the code from the respective tab at the bottom of `agent.py`:

{' ' * (self.list_depth - 1)}- Run synchronously
{' ' * (self.list_depth - 1)}- Stream updates
```
from langchain.messages import HumanMessage

if __name__ == "__main__":
 result = agent.invoke(
 {
 "messages": [
 HumanMessage(
 content="What are the main differences between RAG and fine-tuning for LLM applications?"
 )
 ]
 }
 )

 for msg in result.get("messages", []):
 if hasattr(msg, "content") and msg.content:
 print(msg.content)

```

```
from langchain.messages import HumanMessage
from langgraph.types import Overwrite

if __name__ == "__main__":
 for chunk in agent.stream(
 {
 "messages": [
 HumanMessage(content="Compare Python vs JavaScript for web development")
 ]
 },
 stream_mode="updates",
 ):
 for node, update in chunk.items():
 if not update or not (messages := update.get("messages")):
 continue
 msg_list = messages.value if isinstance(messages, Overwrite) else messages
 for msg in msg_list:
 if hasattr(msg, "content") and msg.content:
 print(msg.content)

```

Run the agent from the project root:

```
python agent.py

```

If you set the `LANGSMITH_API_KEY` environment variable before running, you can view the agent’s traces in [LangSmith](/langsmith/home) to debug and monitor multi-step behavior.

## [​](#full-code)Full code

View the complete [Deep Research example](https://github.com/langchain-ai/deepagents/tree/main/examples/deep_research) on GitHub.

## [​](#next-steps)Next steps

Now that you’ve built the agent, customize it by changing the prompt constants in your agent file to adjust the workflow, delegation strategy, or researcher behavior.
You can also tune the delegation limits to allow for more parallel sub-agents or delegation rounds.
For more information on the concepts in this tutorial, check out the following resources:

{' ' * (self.list_depth - 1)}- [Subagents](/oss/python/deepagents/subagents): Learn how to configure subagents with different tools and prompts

{' ' * (self.list_depth - 1)}- [Customization](/oss/python/deepagents/customization): Customize models, tools, system prompts, and planning behavior

{' ' * (self.list_depth - 1)}- [LangSmith](/langsmith/home): Trace research runs and debug multi-step behavior

{' ' * (self.list_depth - 1)}- [Deep Research Course](https://academy.langchain.com/courses/deep-research-with-langgraph): Full course on deep research with LangGraph

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/deep-research.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a data analysis agentPrevious](/oss/python/deepagents/data-analysis)[Build a content builder agentNext](/oss/python/deepagents/content-builder)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
