# Build a data analysis agent from scratch

[Tutorials](/oss/python/deepagents/data-analysis)[LangChain](/oss/python/langchain/deep-agent-from-scratch)

# Build a data analysis agent from scratch
Copy page

Build a data analysis agent step by step using create_agent and deepagents middleware.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This guide builds a data analysis agent from first principles using `create_agent` and deepagents middleware. Rather than starting with `create_deep_agent`, we assemble the harness one piece at a time: so you can see exactly what each component adds and swap in only what your use case needs.
The agent we’ll build:

{' ' * (self.list_depth - 1)}- Accepts a CSV file for analysis

{' ' * (self.list_depth - 1)}- Writes and executes Python code in an isolated sandbox

{' ' * (self.list_depth - 1)}- Delegates visualization work to a specialized subagent

{' ' * (self.list_depth - 1)}- Loads data analysis patterns from a skills file

## [​](#setup)Setup

```
pip install deepagents langsmith

```

Enable LangSmith tracing to inspect every step:

```
export LANGSMITH_TRACING=true
export LANGSMITH_API_KEY=...

```

---

## [​](#step-1-the-minimal-agent)Step 1: The minimal agent

A model, a loop. Nothing else yet.

```
from langchain.agents import create_agent

agent = create_agent("anthropic:claude-sonnet-4-6", tools=[])

```

This runs, but the agent has no filesystem and no way to execute code. The next steps add those.

---

## [​](#step-2-add-a-sandbox-backend)Step 2: Add a sandbox backend

`LangSmithSandbox` gives the agent an isolated environment with a filesystem and an `execute` tool for running shell commands. The agent can install packages, write scripts, and run them: without touching the host.

```
from langchain.agents import create_agent
from langsmith.sandbox import SandboxClient
from deepagents.backends.langsmith import LangSmithSandbox
from deepagents.middleware import FilesystemMiddleware

client = SandboxClient()
sandbox = client.create_sandbox(template_name="deepagents-deploy")
backend = LangSmithSandbox(sandbox=sandbox)

agent = create_agent(
 "anthropic:claude-sonnet-4-6",
 tools=[],
 middleware=[FilesystemMiddleware(backend=backend)],
)

```

[`FilesystemMiddleware`](https://reference.langchain.com/python/deepagents/middleware/filesystem/FilesystemMiddleware) adds `read_file`, `write_file`, `edit_file`, `glob`, and `grep`. Because `LangSmithSandbox` implements the sandbox protocol, it also adds `execute`: the agent can now run shell commands.
Upload a CSV and invoke:

```
import csv, io

rows = [
 ["Date", "Product", "Units", "Revenue"],
 ["2025-08-01", "Widget A", 10, 250],
 ["2025-08-02", "Widget B", 5, 125],
 ["2025-08-03", "Widget A", 7, 175],
 ["2025-08-04", "Widget C", 3, 90],
]
buf = io.StringIO()
csv.writer(buf).writerows(rows)
backend.upload("sales.csv", buf.getvalue().encode())

result = agent.invoke({
 "messages": [{"role": "user", "content": "Analyze sales.csv. Summarize trends."}]
})

```

---

## [​](#step-3-add-context-management)Step 3: Add context management

For longer analysis sessions the context window fills. `SummarizationMiddleware` compresses history automatically so the agent keeps working without hitting token limits.

```
from deepagents.middleware import FilesystemMiddleware, SummarizationMiddleware

model = "anthropic:claude-sonnet-4-6"

agent = create_agent(
 model=model,
 tools=[],
 middleware=[
 FilesystemMiddleware(backend=backend),
 SummarizationMiddleware(model=model, backend=backend),
 ],
)

```

---

## [​](#step-4-add-skills)Step 4: Add skills

Skills give the agent on-demand domain knowledge via progressive disclosure: loaded only when the current task calls for it. Create a skill file in your skills directory:

```
skills/
 pandas-patterns/
 SKILL.md

```

```
---
name: pandas-patterns
description: Common pandas and matplotlib patterns for data analysis and visualization
---

## Data loading
Use `pd.read_csv()` for CSV files. Always check `df.info()` and `df.describe()` first.

## Visualization
Use `matplotlib` for bar charts, `seaborn` for statistical plots.
Save figures with `plt.savefig("output.png", dpi=150, bbox_inches="tight")`.

## Reporting
Write a markdown summary to `report.md` alongside any generated charts.

```

```
from deepagents.middleware import FilesystemMiddleware, SkillsMiddleware, SummarizationMiddleware

agent = create_agent(
 model=model,
 tools=[],
 middleware=[
 FilesystemMiddleware(backend=backend),
 SummarizationMiddleware(model=model, backend=backend),
 SkillsMiddleware(backend=backend, sources=["./skills/"]),
 ],
)

```

---

## [​](#step-5-add-a-visualization-subagent)Step 5: Add a visualization subagent

Some tasks benefit from isolation. A visualization subagent runs in its own context window, keeping chart generation separate from the main analysis: and enabling parallel execution.

```
from langchain.agents.middleware import TodoListMiddleware
from deepagents import SubAgent
from deepagents.middleware import (
 FilesystemMiddleware,
 SkillsMiddleware,
 SubAgentMiddleware,
 SummarizationMiddleware,
)

visualizer: SubAgent = {
 "name": "visualizer",
 "description": "Generates charts and visualizations from data files in the sandbox.",
 "system_prompt": "You are a data visualization specialist. Write Python scripts using matplotlib and seaborn. Save all figures as PNG files.",
 "tools": [],
}

agent = create_agent(
 model=model,
 tools=[],
 middleware=[
 FilesystemMiddleware(backend=backend),
 SummarizationMiddleware(model=model, backend=backend),
 SkillsMiddleware(backend=backend, sources=["./skills/"]),
 TodoListMiddleware(),
 SubAgentMiddleware(backend=backend, subagents=[visualizer]),
 ],
)

```

The main agent handles analysis and planning; it delegates chart generation to the `visualizer` subagent via the `task` tool.

---

## [​](#what-you-built)What you built

| 
 | Middleware | What it adds
 | [`FilesystemMiddleware`](https://reference.langchain.com/python/deepagents/middleware/filesystem/FilesystemMiddleware) + `LangSmithSandbox` | Isolated filesystem + `execute` tool
 | [`SummarizationMiddleware`](https://reference.langchain.com/python/langchain/agents/middleware/summarization/SummarizationMiddleware) | Automatic context compression
 | [`SkillsMiddleware`](https://reference.langchain.com/python/deepagents/middleware/skills/SkillsMiddleware) | Domain knowledge loaded on demand
 | [`TodoListMiddleware`](https://reference.langchain.com/python/langchain/agents/middleware/todo/TodoListMiddleware) + [`SubAgentMiddleware`](https://reference.langchain.com/python/deepagents/middleware/subagents/SubAgentMiddleware) | Parallel visualization subagent
This is the same foundation as `create_deep_agent`: assembled manually so you control exactly what’s included. The possibilities don’t end here: see [Prebuilt middleware](/oss/python/langchain/middleware/built-in) for the full list of composable capabilities, and the [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) reference for all configuration options.
For a pre-assembled version, see [`create_deep_agent`](https://reference.langchain.com/python/deepagents/graph/create_deep_agent) and [Customize Deep Agents](/oss/python/deepagents/customization). For the full data analysis example using `create_deep_agent`, see [Data analysis](/oss/python/deepagents/data-analysis).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/deep-agent-from-scratch.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a content builder agentPrevious](/oss/python/deepagents/content-builder)[Build a semantic search engine with LangChainNext](/oss/python/langchain/knowledge-base)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
