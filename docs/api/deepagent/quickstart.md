# Quickstart

[Get started](/oss/python/deepagents/quickstart)

# Quickstart
Copy page

Build your first deep agent in minutesCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This guide walks you through creating your first deep agent with planning, file system tools, and subagent capabilities. You’ll build a research agent that can conduct research and write reports.
**Using an AI coding assistant?**

{' ' * (self.list_depth - 1)}- Install the [LangChain Docs MCP server](/use-these-docs) to give your agent access to up-to-date LangChain documentation and examples.

{' ' * (self.list_depth - 1)}- Install [LangChain Skills](https://github.com/langchain-ai/langchain-skills) to improve your agent’s performance on LangChain ecosystem tasks.

## [​](#prerequisites)Prerequisites

Before you begin, make sure you have an API key from a model provider (e.g., Gemini, Anthropic, OpenAI).
Deep Agents require a model that supports [tool calling](/oss/python/langchain/models#tool-calling). See [customization](/oss/python/deepagents/customization#model) for how to configure your model.

## [​](#step-1-install-dependencies)Step 1: Install dependencies

pipuv
```
pip install deepagents tavily-python

```

This guide uses [Tavily](https://tavily.com/) as an example search provider, but you can substitute any search API (e.g., DuckDuckGo, SerpAPI, Brave Search).

## [​](#step-2-set-up-your-api-keys)Step 2: Set up your API keys

{' ' * (self.list_depth - 1)}- Google
{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- OpenRouter
{' ' * (self.list_depth - 1)}- Fireworks
{' ' * (self.list_depth - 1)}- Baseten
{' ' * (self.list_depth - 1)}- Ollama
{' ' * (self.list_depth - 1)}- Other
```
export GOOGLE_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
export OPENAI_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
export ANTHROPIC_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
export OPENROUTER_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
export FIREWORKS_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
export BASETEN_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
# Local: Ollama must be running on your machine
# Cloud: Set your Ollama API key for hosted inference
export OLLAMA_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```

```
# Set the API key for your provider
export <PROVIDER>_API_KEY="your-api-key"
export TAVILY_API_KEY="your-tavily-api-key"

```
Deep Agents work with any [LangChain chat model](/oss/python/deepagents/models#supported-models). Set the API key for your provider.

## [​](#step-3-create-a-search-tool)Step 3: Create a search tool

```
import os
from typing import Literal

from tavily import TavilyClient
from deepagents import create_deep_agent

tavily_client = TavilyClient(api_key=os.environ["TAVILY_API_KEY"])

def internet_search(
 query: str,
 max_results: int = 5,
 topic: Literal["general", "news", "finance"] = "general",
 include_raw_content: bool = False,
):
 """Run a web search"""
 return tavily_client.search(
 query,
 max_results=max_results,
 include_raw_content=include_raw_content,
 topic=topic,
 )

```

## [​](#step-4-create-a-deep-agent)Step 4: Create a deep agent

Pass a `model` string in `provider:model` format, or an [initialized model instance](/oss/python/deepagents/models#configure-model-parameters). See [supported models](/oss/python/deepagents/models#supported-models) for all providers and [suggested models](/oss/python/deepagents/models#suggested-models) for tested recommendations.
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
# System prompt to steer the agent to be an expert researcher
research_instructions = """You are an expert researcher. Your job is to conduct thorough research and then write a polished report.

You have access to an internet search tool as your primary means of gathering information.

## `internet_search`

Use this to run an internet search for a given query. You can specify the max number of results to return, the topic, and whether raw content should be included.
"""

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[internet_search],
 system_prompt=research_instructions,
)

```

## [​](#step-5-run-the-agent)Step 5: Run the agent

```
result = agent.invoke({"messages": [{"role": "user", "content": "What is langgraph?"}]})

# Print the agent's response
print(result["messages"][-1].content)

```

Trace your agent’s planning steps, tool calls, and subagent delegation with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-deepagents-quickstart). Follow the [observability quickstart](/langsmith/observability-quickstart) to get set up.

## [​](#how-does-it-work)How does it work?

Your deep agent automatically:

{' ' * (self.list_depth - 1)}- **Plans its approach** using the built-in [`write_todos`](/oss/python/deepagents/harness#task-planning) tool to break down the research task.

{' ' * (self.list_depth - 1)}- **Conducts research** by calling the `internet_search` tool to gather information.

{' ' * (self.list_depth - 1)}- **Manages context** by using file system tools ([`write_file`](/oss/python/deepagents/harness#virtual-filesystem-access), [`read_file`](/oss/python/deepagents/harness#virtual-filesystem-access)) to offload large search results.

{' ' * (self.list_depth - 1)}- **Spawns subagents** as needed to delegate complex subtasks to specialized subagents.

{' ' * (self.list_depth - 1)}- **Synthesizes a report** to compile findings into a coherent response.

## [​](#examples)Examples

For agents, patterns, and applications you can build with Deep Agents, see [Examples](https://github.com/langchain-ai/deepagents/tree/main/examples).

## [​](#streaming)Streaming

Deep Agents have built-in [streaming](/oss/python/langchain/event-streaming) for real-time updates from agent execution using LangGraph.
This allows you to observe output progressively and review and debug agent and subagent work, such as tool calls, tool results, and LLM responses.

## [​](#next-steps)Next steps

Now that you’ve built your first deep agent:

{' ' * (self.list_depth - 1)}- **Customize your agent**: Learn about [customization options](/oss/python/deepagents/customization), including custom system prompts, tools, and subagents.

{' ' * (self.list_depth - 1)}- **Add long-term memory**: Enable [persistent memory](/oss/python/deepagents/memory) across conversations.

{' ' * (self.list_depth - 1)}- **Deploy to production**: Use [Managed Deep Agents](/langsmith/deploy-managed-deep-agent) to create, run, and operate deep agents in LangSmith.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/quickstart.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Deep Agents overviewPrevious](/oss/python/deepagents/overview)[Customize Deep AgentsNext](/oss/python/deepagents/customization)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
