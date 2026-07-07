# Skills

[Core capabilities](/oss/python/deepagents/harness)

# Skills
Copy page

Learn how to extend your deep agent’s capabilities with skillsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Skills are reusable agent capabilities that provide specialized workflows and domain knowledge.
You can use [Agent Skills](https://agentskills.io/) to provide your deep agent with new capabilities and expertise. For ready-to-use skills that improve your agent’s performance on LangChain ecosystem tasks, see the [LangChain Skills](https://github.com/langchain-ai/langchain-skills) repository.
Deep agent skills follow the [Agent Skills specification](https://agentskills.io/specification) and add additional capability for interpreter skills, which makes it possible to provide skills with importable functions that an interpreter can call.

## [​](#what-are-skills)What are skills

Skills are a directory of folders, where each folder has one or more files that contain context the agent can use:

{' ' * (self.list_depth - 1)}- A `SKILL.md` file containing instructions and metadata about the skill

{' ' * (self.list_depth - 1)}- Additional scripts (optional)

{' ' * (self.list_depth - 1)}- Additional reference info, such as docs (optional)

{' ' * (self.list_depth - 1)}- Additional assets, such as templates and other resources (optional)

Any additional assets (scripts, docs, templates, or other resources) must be referenced in the `SKILL.md` file with information on what the file contains and how to use it so the agent can decide when to use them.

## [​](#how-skills-work)How skills work

When you create a deep agent, you can pass in a list of directories containing skills. As the agent starts, it reads through the frontmatter of each `SKILL.md` file.
When the agent receives a prompt, the agent checks if it can use any skills while fulfilling the prompt. If it finds a matching prompt, it then reviews the rest of the skill files. This pattern of only reviewing the skill information when needed is called *progressive disclosure*.

## [​](#example)Example

You might have a skills folder that contains a skill to use a docs site in a certain way, as well as another skill to search the arXiv preprint repository of research papers:

```
 skills/
 ├── langgraph-docs
 │ └── SKILL.md
 └── arxiv_search
 ├── SKILL.md
 └── arxiv_search.py # code for searching arXiv

```

The `SKILL.md` file always follows the same pattern, starting with metadata in the frontmatter and followed by the instructions for the skill.
The following example shows a skill that gives instructions on how to provide relevant langgraph docs when prompted:

```
---
name: langgraph-docs
description: Use this skill for requests related to LangGraph in order to fetch relevant documentation to provide accurate, up-to-date guidance.
module: index.ts
---

# langgraph-docs

## Overview

This skill explains how to access LangGraph Python documentation to help answer questions and guide implementation.

## Instructions

### 1. Fetch the Documentation Index

Use the fetch_url tool to read the following URL:
https://docs.langchain.com/llms.txt

This provides a structured list of all available documentation with descriptions.

### 2. Select Relevant Documentation

Based on the question, identify 2-4 most relevant documentation URLs from the index. Prioritize:

- Specific how-to guides for implementation questions
- Core concept pages for understanding questions
- Tutorials for end-to-end examples
- Reference docs for API details

### 3. Fetch Selected Documentation

Use the fetch_url tool to read the selected documentation URLs.

### 4. Provide accurate guidance

After reading the documentation, answer the user's question using the relevant LangGraph docs you fetched.

In your response:

- Give a direct answer first.
- Include the minimum necessary context and any key steps or API names.
- Avoid quoting long passages. Paraphrase and link instead.

### 5. Provide the regular links for the used references

At the end of your response, include a **References** section listing the page URLs you used.

`llms.txt` uses Markdown link targets that typically end in `.md`. Use the helper from this skill module to resolve those into the actual page URLs before listing them as references.

```typescript
const { resolveLlmsUrl } = await import("@/skills/langgraph-docs");

// llms.txt uses Markdown link targets that typically end in `.md`.
// Convert those into the actual page URLs before fetching.
const llmsUrls = [
 "https://docs.langchain.com/oss/langgraph/concepts.md",
 "https://docs.langchain.com/oss/langgraph/concepts.md",
 "https://docs.langchain.com/oss/langgraph/tutorials.md",
];

const pageUrls = [...new Set(llmsUrls.map(resolveLlmsUrl))];
pageUrls;
```

```

The referenced helper code would be placed in `index.ts`:

```
// index.ts
export function resolveLlmsUrl(url: string) {
 return url.endsWith(".md") ? url.slice(0, -3) : url;
}

```

For more example skills, see [Deep Agents example skills](https://github.com/langchain-ai/deepagents/tree/main/libs/cli/examples/skills).
**Important**Refer to the full [Agent Skills Specification](https://agentskills.io/specification) for information on constraints and best practices when authoring skill files. Notably:

{' ' * (self.list_depth - 1)}- The `description` field is truncated to 1024 characters if it exceeds that length.

{' ' * (self.list_depth - 1)}- In Deep Agents, `SKILL.md` files must be under 10 MB. Files exceeding this limit are skipped during skill loading.

### [​](#full-example)Full example

The following example shows a `SKILL.md` file using all available frontmatter fields:

```
---
name: langgraph-docs
description: Use this skill for requests related to LangGraph in order to fetch relevant documentation to provide accurate, up-to-date guidance.
license: MIT
compatibility: Requires internet access for fetching documentation URLs
metadata:
 author: langchain
 version: "1.0"
allowed-tools: fetch_url
module: index.ts
---

# langgraph-docs

## Overview

This skill explains how to access LangGraph Python documentation to help answer questions and guide implementation.

## Instructions

### 1. Fetch the documentation index

Use the fetch_url tool to read the following URL:
https://docs.langchain.com/llms.txt

This provides a structured list of all available documentation with descriptions.

### 2. Select relevant documentation

Based on the question, identify 2-4 most relevant documentation URLs from the index. Prioritize:

- Specific how-to guides for implementation questions
- Core concept pages for understanding questions
- Tutorials for end-to-end examples
- Reference docs for API details

### 3. Fetch selected documentation

Use the fetch_url tool to read the selected documentation URLs.

### 4. Provide accurate guidance

After reading the documentation, answer the user's question using the relevant LangGraph docs you fetched.

In your response:

- Give a direct answer first.
- Include the minimum necessary context and any key steps or API names.
- Avoid quoting long passages. Paraphrase and link instead.

### 5. Provide the regular links for the used references

At the end of your response, include a **References** section listing the page URLs you used.

`llms.txt` uses Markdown link targets that typically end in `.md`. Use the helper from this skill module to resolve those into the actual page URLs before listing them as references.

```typescript
const { resolveLlmsUrl } = await import("@/skills/langgraph-docs");

// llms.txt uses Markdown link targets that typically end in `.md`.
// Convert those into the actual page URLs before fetching.
const llmsUrls = [
 "https://docs.langchain.com/oss/langgraph/concepts.md",
 "https://docs.langchain.com/oss/langgraph/concepts.md",
 "https://docs.langchain.com/oss/langgraph/tutorials.md",
];

const pageUrls = [...new Set(llmsUrls.map(resolveLlmsUrl))];
pageUrls;
```

```
See all 70 lines

## [​](#usage)Usage

Pass the skills directory when creating your deep agent:

{' ' * (self.list_depth - 1)}- StateBackend
{' ' * (self.list_depth - 1)}- StoreBackend
{' ' * (self.list_depth - 1)}- FilesystemBackendGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from urllib.request import urlopen
from deepagents import create_deep_agent
from deepagents.backends import StateBackend
from deepagents.backends.utils import create_file_data
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.checkpoint.memory import MemorySaver

checkpointer = MemorySaver()
backend = StateBackend()

skill_url = "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/libs/cli/examples/skills/langgraph-docs/SKILL.md"
with urlopen(skill_url) as response:
 skill_content = response.read().decode('utf-8')

skills_files = {
 "/skills/langgraph-docs/SKILL.md": create_file_data(skill_content),
}

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 skills=["/skills/"],
 checkpointer=checkpointer,
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)], # for interpreter skills
)

result = agent.invoke(
 {
 "messages": [{"role": "user", "content": "What is langgraph?"}],
 # Seed the default StateBackend's in-state filesystem (virtual paths must start with "/").
 "files": skills_files,
 },
 config={"configurable": {"thread_id": "12345"}},
)

```

```
from urllib.request import urlopen
from deepagents import create_deep_agent
from deepagents.backends import StoreBackend
from deepagents.backends.utils import create_file_data
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.store.memory import InMemoryStore

store = InMemoryStore()
backend = StoreBackend(namespace=lambda _rt: ("filesystem",))

skill_url = "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/libs/cli/examples/skills/langgraph-docs/SKILL.md"
with urlopen(skill_url) as response:
 skill_content = response.read().decode('utf-8')

store.put(
 namespace=("filesystem",),
 key="/skills/langgraph-docs/SKILL.md",
 value=create_file_data(skill_content),
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 store=store,
 skills=["/skills/"],
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)],
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What is langgraph?"}]},
 config={"configurable": {"thread_id": "12345"}},
)

```

```
from deepagents import create_deep_agent
from deepagents.backends.filesystem import FilesystemBackend
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.checkpoint.memory import MemorySaver

# Checkpointer is REQUIRED for human-in-the-loop
checkpointer = MemorySaver()
root_dir = "/Users/user/{project}"
backend = FilesystemBackend(root_dir=root_dir)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 skills=[str(Path(root_dir) / "skills")],
 interrupt_on={
 "write_file": True,
 "read_file": False,
 "edit_file": True,
 },
 checkpointer=checkpointer, # Required!
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)], # for interpreter skills
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What is langgraph?"}]},
 config={"configurable": {"thread_id": "12345"}},
)

```

[​](#param-skills)skillslist[str]List of skill source paths.Paths must be specified using forward slashes and are relative to the backend’s root.

{' ' * (self.list_depth - 1)}- If omitted, no skills are loaded.

{' ' * (self.list_depth - 1)}- When using `StateBackend` (default), provide skill files with `invoke(files={...})`. Use `create_file_data()` from `deepagents.backends.utils` to format file contents; raw strings are not supported.

{' ' * (self.list_depth - 1)}- With `FilesystemBackend`, skills are loaded from disk relative to the backend’s `root_dir`.
Later sources override earlier ones for skills with the same name (last one wins).
The SDK only loads the sources you pass in `skills`. It does not automatically scan CLI directories such as `~/.deepagents/...` or `~/.agents/...`.For CLI storage conventions, see [App data](/oss/python/deepagents/data-locations).

Emulating CLI source order in SDKIf you want CLI-style layering in SDK code, pass all desired sources explicitly in lowest-to-highest precedence order:
```
[
"<user-home>/.deepagents/{agent}/skills/",
"<user-home>/.agents/skills/",
"<project-root>/.deepagents/skills/",
"<project-root>/.agents/skills/",
]

```
Then pass that ordered list as `skills` when creating your agent.

## [​](#source-precedence)Source precedence

When multiple skill sources contain a skill with the same name, the skill from the source listed later in the `skills` array takes precedence (last one wins). This lets you layer skills from different origins.

```
# If both sources contain a skill named "web-search",
# the one from "/skills/project/" wins (loaded last).
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 skills=["/skills/user/", "/skills/project/"],
 ...
)

```

## [​](#skills-for-subagents)Skills for subagents

When you use [subagents](/oss/python/deepagents/subagents), you can configure which skills each type has access to:

{' ' * (self.list_depth - 1)}- **General-purpose subagent**: Automatically inherits skills from the main agent when you pass `skills` to `create_deep_agent`. No additional configuration is needed.

{' ' * (self.list_depth - 1)}- **Custom subagents**: Do not inherit the main agent’s skills. Add a `skills` parameter to each subagent definition with that subagent’s skill source paths.

Skill state is fully isolated: the main agent’s skills are not visible to subagents, and subagent skills are not visible to the main agent.

```
from deepagents import create_deep_agent

research_subagent = {
 "name": "researcher",
 "description": "Research assistant with specialized skills",
 "system_prompt": "You are a researcher.",
 "tools": [web_search],
 "skills": ["/skills/research/", "/skills/web-search/"], # Subagent-specific skills
}

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 skills=["/skills/main/"], # Main agent and GP subagent get these
 subagents=[research_subagent], # Researcher gets only its own skills
)

```

For more information on subagent configuration and skills inheritance, see [Subagents](/oss/python/deepagents/subagents).

## [​](#what-the-agent-sees)What the agent sees

When skills are configured, a “Skills System” section is injected into the agent’s system prompt. The agent uses this information to follow a three-step process:

{' ' * (self.list_depth - 1)}- **Match**—When a user prompt arrives, the agent checks whether any skill’s description matches the task.

{' ' * (self.list_depth - 1)}- **Read**—If a skill applies, the agent reads the full `SKILL.md` file using the path shown in its skills list.

{' ' * (self.list_depth - 1)}- **Execute**—The agent follows the skill’s instructions and accesses any supporting files (scripts, templates, reference docs) as needed.

Write clear, specific descriptions in your `SKILL.md` frontmatter. The agent decides whether to use a skill based on the description alone—detailed descriptions lead to better skill matching.

## [​](#execute-code-with-skills)Execute code with skills

Skills support code execution in two ways:

{' ' * (self.list_depth - 1)}- [Use interpreter skills](#use-interpreter-skills) when the agent needs reusable, importable helpers available from interpreter code.

{' ' * (self.list_depth - 1)}- [Execute skill scripts in a sandbox](#execute-skill-scripts-in-a-sandbox) when the agent needs to install dependencies, run tests, call CLIs, or work with an operating-system filesystem.

### [​](#use-interpreter-skills)Use interpreter skills

Interpreter skills are skills that expose code modules to an [interpreter](/oss/python/deepagents/interpreters). Regular skills give the agent instructions and context. Interpreter skills also give the agent importable functions it can call from interpreter code.
This lets you package domain-specific logic once and make it available as a deterministic building block inside the agent’s workspace. Instead of asking the model to re-create a parser, scorer, normalizer, validator, or aggregation routine from scratch, the agent can import a tested helper and compose it with tools, subagents, and runtime state.
Use interpreter skills for code that should be:

{' ' * (self.list_depth - 1)}- **Reusable** across prompts, agents, or projects.

{' ' * (self.list_depth - 1)}- **Deterministic** enough that you want the same behavior every time.

{' ' * (self.list_depth - 1)}- **Too detailed** to keep in the model context as instructions.

{' ' * (self.list_depth - 1)}- **Useful inside larger workflows**, such as scoring search results, normalizing API responses, validating records, grouping rows, or converting data into a report-ready shape.

To make a skill importable:
1[](#)

#### Add a module entry
Add a `module` key to the skill’s `SKILL.md` frontmatter. The value is a JavaScript or TypeScript file path relative to the skill directory.2[](#)

#### Configure skills normally
Pass the skill source path with the `skills` argument when creating the agent.3[](#)

#### Use the same backend
Configure the interpreter middleware with the same backend that `SkillsMiddleware` uses to load skill files.4[](#)

#### Import from interpreter code
The agent imports the helper module with `await import("@/skills/<name>")`.
Minimal skill layout:

```
skills/
`-- order-helpers/
 |-- SKILL.md
 `-- index.ts

```

```
---
name: order-helpers
description: Helper functions for normalizing and grouping order records.
module: index.ts
---

# order-helpers

Use this skill when order records need deterministic cleanup or aggregation.

Import these utilities into the REPL in order to interact with order data:

```typescript
const { groupByStatus } = await import("@/skills/order-helpers");
groupByStatus(...);
```

```

```
// skills/order-helpers/index.ts
interface Order {
 id: string;
 status: string;
}

export function groupByStatus(orders: Order[]) {
 return orders.reduce((acc, order) => {
 acc[order.status] = acc[order.status] ?? [];
 acc[order.status].push(order);
 return acc;
 }, {});
}

```

Then configure the agent:

```
from deepagents import create_deep_agent
from deepagents.backends import StateBackend
from langchain_quickjs import CodeInterpreterMiddleware

backend = StateBackend()

agent = create_deep_agent(
 model="openai:gpt-5.4",
 backend=backend,
 skills=["/skills/"],
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)],
)

```

The agent can now import the module from interpreter code:

```
const { groupByStatus } = await import("@/skills/order-helpers");

const grouped = groupByStatus(orders);
grouped;

```

### [​](#execute-skill-scripts-in-a-sandbox)Execute skill scripts in a sandbox

Skills can include scripts alongside the `SKILL.md` file, such as, for example, a Python file that performs a search or data transformation. The agent can *read* these scripts from any backend, but to *execute* them, the agent needs access to a shell — which only [sandbox backends](/oss/python/deepagents/sandboxes) provide.
When you use a [CompositeBackend](https://reference.langchain.com/python/deepagents/backends/composite/CompositeBackend) that routes skills to a [StoreBackend](https://reference.langchain.com/python/deepagents/backends/store/StoreBackend) for persistence while using a sandbox as the default backend, skill files live in the store rather than in the sandbox is where code runs. For sandboxes to be able to use the scripts, you must use [custom middleware](/oss/python/langchain/middleware/custom) to upload skill scripts into the sandbox before the agent starts:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
import asyncio
from pathlib import Path
from typing import Any

from daytona import Daytona
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StoreBackend
from deepagents.backends.utils import create_file_data
from langchain.agents.middleware import AgentMiddleware, AgentState

from langchain_daytona import DaytonaSandbox
from langgraph.runtime import Runtime
from langgraph.store.memory import InMemoryStore

# Identical skill bundles for every user: one shared store namespace.
SKILLS_SHARED_NAMESPACE = ("skills", "builtin")

class SkillSandboxSyncMiddleware(AgentMiddleware[AgentState, Any, Any]):
 """Copy shared skill files from the store into the sandbox before each agent run."""

 def __init__(self, backend: CompositeBackend) -> None:
 super().__init__()
 self.backend = backend

 async def abefore_agent(self, state: AgentState, runtime: Runtime[Any]) -> None:
 store = runtime.store

 files: list[tuple[str, bytes]] = []
 for item in await store.asearch(SKILLS_SHARED_NAMESPACE):
 key = str(item.key)
 if ".." in key or any(c in key for c in ("*", "?")):
 msg = f"Invalid key: {key}"
 raise ValueError(msg)
 normalized = key if key.startswith("/") else f"/{key}"
 # CompositeBackend routes paths and batches uploads to the right backend.
 files.append((f"/skills{normalized}", item.value["content"].encode()))

 if files:
 await self.backend.aupload_files(files)

async def seed_skill_store(store: InMemoryStore) -> None:
 """Load canonical skill files from disk into the shared store namespace (run once at deploy).
 You can retrieve skills from any source (local filesystem, remote URL, etc.).
 """
 skills_dir = Path(__file__).resolve().parent / "skills"
 for file_path in sorted(p for p in skills_dir.rglob("*") if p.is_file()):
 rel = file_path.relative_to(skills_dir).as_posix()
 key = f"/{rel}"
 await store.aput(
 SKILLS_SHARED_NAMESPACE,
 key,
 create_file_data(file_path.read_text(encoding="utf-8")),
 )

async def main() -> None:
 store = InMemoryStore()
 await seed_skill_store(store)

 daytona = Daytona()
 sandbox = daytona.create()
 sandbox_backend = DaytonaSandbox(sandbox=sandbox)

 backend = CompositeBackend(
 default=sandbox_backend,
 routes={
 "/skills/": StoreBackend(
 store=store,
 namespace=lambda _rt: SKILLS_SHARED_NAMESPACE,
 ),
 },
 )

 try:
 agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 skills=["/skills/"],
 store=store,
 middleware=[SkillSandboxSyncMiddleware(backend)],
 )

 finally:
 sandbox.stop()

if __name__ == "__main__":
 asyncio.run(main())

```

The middleware’s `before_agent` hook runs before each agent invocation, reading skill files from that shared namespace and uploading them into the sandbox filesystem. Once synced, the agent can execute scripts with the `execute` tool just like any other file in the sandbox.
For a more complete example that also syncs [memories](/oss/python/deepagents/memory) bidirectionally, see [syncing skills and memories with custom middleware](/oss/python/deepagents/going-to-production#example-syncing-skills-and-memories-with-custom-middleware).

## [​](#skills-vs-memory)Skills vs. memory

Skills and [memory](/oss/python/deepagents/customization#memory) (`AGENTS.md` files) serve different purposes:

| 
 | | Skills | Memory
 | **Purpose** | On-demand capabilities discovered through progressive disclosure | Persistent context always loaded at startup
 | **Loading** | Read only when the agent determines relevance | Always injected into system prompt
 | **Format** | `SKILL.md` in named directories | `AGENTS.md` files
 | **Layering** | User → project (last wins) | User → project (combined)
 | **Use when** | Instructions are task-specific and potentially large | Context is always relevant (project conventions, preferences)

## [​](#when-to-use-skills-and-tools)When to use skills and tools

These are a few general guidelines for using tools and skills:

{' ' * (self.list_depth - 1)}- Use skills when there is a lot of context to reduce the number of tokens in the system prompt.

{' ' * (self.list_depth - 1)}- Use skills to bundle capabilities together into larger actions and provide additional context beyond single tool descriptions.

{' ' * (self.list_depth - 1)}- Use tools if the agent does not have access to the file system.

Trace how your agent discovers and executes skills with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-deepagents-skills). Follow the [observability quickstart](/langsmith/observability-quickstart) to get set up.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/skills.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[MemoryPrevious](/oss/python/deepagents/memory)[SandboxesNext](/oss/python/deepagents/sandboxes)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
