# Deep Agents overview

# Deep Agents overview
Copy page

Build agents that can plan, use subagents, and leverage file systems for complex tasksCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.The easiest way to start building agents and applications powered by LLMs—with built-in capabilities for task planning, file systems for context management, subagent-spawning, and long-term memory.
You can use deep agents for any task, including complex, multi-step tasks.
Deep Agents is an [“agent harness”](/oss/python/concepts/products#agent-harnesses-like-the-deep-agents-sdk). It is the same core tool calling loop as other agent frameworks, but with built-in capabilities that make agents reliable for real tasks:

## Take actions in an environment
Take actions via tools, read and write files, execute code

## Connect to your data
Load memories, skills, and domain knowledge at the right moment

## Manage growing context
Summarize history and offload large results across long runs

## Parallelize tasks
Delegate to general or specialized subagents running in isolated context windows

## Stay in the loop
Pause for human approval at critical decision points

## Improve over time
Update memory, skills, and prompts based on real usage
See [Harness capabilities](/oss/python/deepagents/harness) for a full breakdown of each component.
The [LangSmith Engine](/langsmith/engine) detects issues in your Deep Agents traces and proposes fixes. You can open a pull request with the proposed fix directly from the Engine tab.
[`deepagents`](https://pypi.org/project/deepagents/) is a standalone library built on top of [LangChain](/oss/python/langchain)’s core building blocks for agents. It uses the [LangGraph](/oss/python/langgraph) runtime for durable execution, streaming, human-in-the-loop, and other features.
The [`deepagents` repository](https://github.com/langchain-ai/deepagents) contains:

{' ' * (self.list_depth - 1)}- **Deep Agents SDK**: A package for building agents that can handle any task

{' ' * (self.list_depth - 1)}- [**Deep Agents Code**](/oss/python/deepagents/code): A terminal coding agent built on the Deep Agents SDK

{' ' * (self.list_depth - 1)}- [**ACP integration**](/oss/python/deepagents/acp): An Agent Client Protocol connector for using deep agents in code editors like Zed

[LangChain](/oss/python/langchain) is the framework that provides the core building blocks for your agents.
To learn more about the differences between LangChain, LangGraph, and Deep Agents, see [Frameworks, runtimes, and harnesses](/oss/python/concepts/products). For a side-by-side comparison with Anthropic’s harness, see [Deep Agents vs. Claude Agent SDK](/oss/python/deepagents/comparison).

## [​](#create-a-deep-agent) Create a deep agent

{' ' * (self.list_depth - 1)}- Google
{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- OpenRouter
{' ' * (self.list_depth - 1)}- Fireworks
{' ' * (self.list_depth - 1)}- Baseten
{' ' * (self.list_depth - 1)}- Ollama
```
# pip install -qU deepagents langchain-google-genai
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-openai
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="openai:gpt-5.4",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-anthropic
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="anthropic:claude-sonnet-4-6",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-openrouter
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="openrouter:anthropic/claude-sonnet-4-6",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-fireworks
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="fireworks:accounts/fireworks/models/qwen3p5-397b-a17b",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-baseten
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="baseten:zai-org/GLM-5",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

```
# pip install -qU deepagents langchain-ollama
from deepagents import create_deep_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_deep_agent(
 model="ollama:devstral-2",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "what is the weather in sf"}]}
)

```

See the [Quickstart](/oss/python/deepagents/quickstart) and [Customization guide](/oss/python/deepagents/customization) to get started building your own agents and applications with Deep Agents.
Trace requests, debug agent behavior, and evaluate outputs with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-deepagents-overview). Follow the [tracing quickstart](/langsmith/trace-with-langchain) to get set up. When ready for production, [deploy to LangSmith Cloud](/langsmith/deploy-to-cloud) for managed hosting.

## [​](#when-to-use-deep-agents)When to use Deep Agents

Use the **Deep Agents SDK** when you want to build agents that can:

{' ' * (self.list_depth - 1)}- **Handle complex, multi-step tasks** that require planning and decomposition

{' ' * (self.list_depth - 1)}- **Manage large amounts of context** through file system tools and [summarization](/oss/python/deepagents/context-engineering#summarization)

{' ' * (self.list_depth - 1)}- **Swap filesystem backends** to use in-memory state, local disk, durable stores, [sandboxes](/oss/python/deepagents/sandboxes), or [your own custom backend](/oss/python/deepagents/backends)

{' ' * (self.list_depth - 1)}- **Execute shell commands** via the `execute` tool when using a [sandbox backend](/oss/python/deepagents/sandboxes)

{' ' * (self.list_depth - 1)}- **Run interpreter code** with [interpreters](/oss/python/deepagents/interpreters) for tool composition, subagent orchestration, and structured data transformations

{' ' * (self.list_depth - 1)}- **Delegate work** to specialized subagents for context isolation

{' ' * (self.list_depth - 1)}- **Persist memory** across conversations and threads

{' ' * (self.list_depth - 1)}- **Control filesystem access** with declarative [permission rules](/oss/python/deepagents/permissions) that restrict which files agents can read or write

{' ' * (self.list_depth - 1)}- **Require human approval** for sensitive operations with [human-in-the-loop](/oss/python/deepagents/human-in-the-loop) workflows

{' ' * (self.list_depth - 1)}- **Use any model** through [provider-agnostic model support](/oss/python/deepagents/models)

For building simpler agents, consider using LangChain’s [`create_agent`](/oss/python/langchain/agents) or building a custom [LangGraph](/oss/python/langgraph/overview) workflow.

## [​](#core-capabilities)Core capabilities

Deep Agents gives an agent the scaffolding it needs to work across long-running, multi-step tasks. The capabilities fall into a few categories:

| 
 | Question | Capability | What it gives the agent
 | What should happen next? | Planning | Break a request into steps, track progress, and revise the plan as new information appears.
 | Where does the work go? | Context management | Store files, intermediate outputs, and large tool results outside the model context window.
 | How does work get done? | Tool use, shell execution, and interpreters | Call your tools, run shell commands in sandboxes, and execute small programs for data transformation or orchestration.
 | How does the agent go deeper without losing focus? | Subagents | Delegate focused work to isolated agents and return only the result to the main agent.
 | How does the agent learn over time? | Memory and skills | Persist useful information across threads and load reusable procedures when they match the task.
 | How do you control risk? | Backends, permissions, and human approval | Choose where files live, limit what the agent can read or write, and pause sensitive actions for review.
These pieces work together as a harness around the model:

{' ' * (self.list_depth - 1)}- **Plan work:** A built-in [`write_todos`](/oss/python/langchain/middleware/built-in#to-do-list) tool helps agents decompose complex tasks, track progress, and adapt their plan.

{' ' * (self.list_depth - 1)}- **Manage context:** File system tools ([`ls`](/oss/python/deepagents/harness#virtual-filesystem-access), [`read_file`](/oss/python/deepagents/harness#virtual-filesystem-access), [`write_file`](/oss/python/deepagents/harness#virtual-filesystem-access), [`edit_file`](/oss/python/deepagents/harness#virtual-filesystem-access)) let agents move large results into files. [Summarization](/oss/python/deepagents/context-engineering#summarization) compacts older conversation history during long runs.

{' ' * (self.list_depth - 1)}- **Run code and commands:** [Sandbox backends](/oss/python/deepagents/sandboxes) provide an `execute` tool for tests, builds, git operations, and system tasks. [Interpreters](/oss/python/deepagents/interpreters) run JavaScript in an in-memory runtime for tool composition, subagent orchestration, and structured data transformations.

{' ' * (self.list_depth - 1)}- **Delegate focused work:** A built-in `task` tool lets agents spawn [subagents](/oss/python/deepagents/subagents) for isolated subtasks, keeping the main agent focused on coordination.

{' ' * (self.list_depth - 1)}- **Remember and reuse knowledge:** [Long-term memory](/oss/python/deepagents/memory) persists information across threads, and [skills](/oss/python/deepagents/skills) package reusable workflows, domain knowledge, and instructions.

{' ' * (self.list_depth - 1)}- **Control execution:** [Pluggable backends](/oss/python/deepagents/backends) determine where files are stored. [Permission rules](/oss/python/deepagents/permissions) restrict filesystem access, and [human-in-the-loop](/oss/python/deepagents/human-in-the-loop) workflows require approval for sensitive tool operations.

{' ' * (self.list_depth - 1)}- **Start from useful defaults:** Deep Agents includes system prompts that guide models to plan before acting, verify work, and manage context. You can customize or replace the defaults.

## [​](#get-started)Get started

## Quickstart
Build your first deep agent

## Customization
Learn about customization options

## Models
Configure models and providers

## Backends
Choose and configure pluggable filesystem backends

## Sandboxes
Execute code in isolated environments

## Interpreters
Compose tools and transform data in QuickJS

## Permissions
Control filesystem access with permission rules

## Human-in-the-loop
Configure approval for sensitive operations

## Code
Use Deep Agents Code

## ACP
Use deep agents in code editors via ACP

## Reference
See the `deepagents` API reference

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/overview.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[QuickstartNext](/oss/python/deepagents/quickstart)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
