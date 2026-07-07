# Harness capabilities

[Core capabilities](/oss/python/deepagents/harness)

# Harness capabilities
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.A Deep Agents harness provides four categories of built-in capabilities that make building long-running, reliable agents easier:

## Execution environment
Tools, virtual filesystem, optional sandbox, and REPL (interpreter)

## Context management
Skills, memory, summarization, context offloading, and prompt caching

## Delegation
Subagent spawning and task planning

## Steering
Human-in-the-loop approval and interrupts
Alongside these four components, [harness profiles](#harness-profiles) let you package per-model configuration into reusable bundles.

![The Deep Agents open harness: planning, virtual filesystem, permissions, subagents, context management, code execution, human-in-the-loop, skills, and memory](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/open-harness.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=aad6b98dc01a1401c96c46e36f3c4dd9)

## [​](#execution-environment)Execution environment

The execution environment is where an agent acts. It has four layers:

{' ' * (self.list_depth - 1)}- **[Tools](#tools)**: custom functions, APIs, and databases the agent can call

{' ' * (self.list_depth - 1)}- **[Virtual filesystem](#virtual-filesystem-access)**: file tools backed by pluggable backends

{' ' * (self.list_depth - 1)}- **[Filesystem permissions](#filesystem-permissions)**: declarative access control over which paths agents can read or write

{' ' * (self.list_depth - 1)}- **[Code execution](#code-execution)**: sandboxed shell execution and an in-process JavaScript interpreter

### [​](#tools)Tools

Pass any Python callable, LangChain tool, or tool dict to `create_deep_agent` via the `tools=` parameter. These are the domain-specific actions your agent can take—web search, database queries, API calls, or any function you define.

```
from deepagents import create_deep_agent

agent = create_deep_agent(
 model="anthropic:claude-sonnet-4-6",
 tools=[search, fetch_page, run_query],
)

```

For more information on defining and using tools, see [Tools](/oss/python/langchain/tools).

### [​](#virtual-filesystem-access)Virtual filesystem access

The harness provides a configurable virtual filesystem which can be backed by different pluggable backends.
The backends support the following file system operations:

| 
 | Tool | Description
 | `ls` | List files in a directory with metadata (size, modified time)
 | `read_file` | Read file contents with line numbers, supports offset/limit for large files. Also supports returning multimodal content blocks for non-text files (images, video, audio, and documents). See supported extensions below.
 | `write_file` | Create new files
 | `edit_file` | Perform exact string replacements in files (with global replace mode)
 | `glob` | Find files matching patterns (e.g., `**/*.py`)
 | `grep` | Search file contents with multiple output modes (files only, content with context, or counts)
 | `execute` | Run shell commands in the environment (available with [sandbox backends](/oss/python/deepagents/sandboxes) only)

Supported multimodal file extensions

| 
 | Type | Extensions
 | [Image](/oss/python/langchain/messages#multimodal) | `.png`, `.jpg`, `.jpeg`, `.gif`, `.webp`, `.heic`, `.heif`
 | [Video](/oss/python/langchain/messages#multimodal) | `.mp4`, `.mpeg`, `.mov`, `.avi`, `.flv`, `.mpg`, `.webm`, `.wmv`, `.3gpp`
 | [Audio](/oss/python/langchain/messages#multimodal) | `.wav`, `.mp3`, `.aiff`, `.aac`, `.ogg`, `.flac`
 | [File](/oss/python/langchain/messages#multimodal) | `.pdf`, `.ppt`, `.pptx`

Running without the default filesystem toolsTo hide the filesystem tools listed above from the model, register a [harness profile](#harness-profiles) with `excluded_tools`:
```
from deepagents import HarnessProfile, register_harness_profile

register_harness_profile(
 "anthropic:claude-sonnet-4-6",
 HarnessProfile(
 excluded_tools=frozenset(
 {"ls", "read_file", "write_file", "edit_file", "glob", "grep"}
 ),
 ),
)

```
Removing `FilesystemMiddleware` itself via `excluded_middleware` is intentionally rejected—use `excluded_tools` to hide only the model-visible tool surface and leave the middleware in place. To remove the `task` tool, see [Running without subagents](/oss/python/deepagents/subagents#running-without-subagents).
The virtual filesystem is used by several other harness capabilities such as skills, memory, code execution, and context management.
You can also use the file system when building custom tools and middleware for Deep Agents.
For more information, see [backends](/oss/python/deepagents/backends).

### [​](#filesystem-permissions)Filesystem permissions

The harness supports declarative permission rules that control which files and directories the agent can read or write. Permissions apply to the built-in filesystem tools listed above and are evaluated in declaration order with first-match-wins semantics.
**How it works:**

{' ' * (self.list_depth - 1)}- Pass a list of rules to `permissions=` when creating the agent

{' ' * (self.list_depth - 1)}- Each rule specifies `operations` (`"read"`, `"write"`), `paths` (glob patterns), and `mode` (`"allow"` or `"deny"`)

{' ' * (self.list_depth - 1)}- The first matching rule wins. If no rule matches, the operation is allowed.

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Restrict agents to specific directories (e.g., `/workspace/`)

{' ' * (self.list_depth - 1)}- Protect sensitive files (e.g., `.env`, credentials)

{' ' * (self.list_depth - 1)}- Give subagents narrower access than the parent agent

Permissions do not apply to [sandbox backends](/oss/python/deepagents/sandboxes), which support arbitrary command execution via the `execute` tool. For custom validation logic, use [backend policy hooks](/oss/python/deepagents/backends#add-policy-hooks).
For the full rule structure, examples, and subagent inheritance, see [Permissions](/oss/python/deepagents/permissions).

### [​](#code-execution)Code execution

Deep Agents supports code execution in two ways:

{' ' * (self.list_depth - 1)}- [Sandbox backends](/oss/python/deepagents/sandboxes) expose an `execute` tool for shell commands in an isolated environment.

{' ' * (self.list_depth - 1)}- [Interpreters](/oss/python/deepagents/interpreters) add an `eval` tool that runs JavaScript in a scoped QuickJS runtime.

Use sandbox backends when the agent needs to install dependencies, run tests, call CLIs, or work with an operating-system filesystem. Sandbox backends implement the `SandboxBackendProtocolV2`; when detected, the harness adds the `execute` tool to the agent’s available tools.
Use interpreters when the agent needs a lightweight programmable layer for loops, batching, deterministic data transformations, or programmatic tool calling. Interpreters do not provide shell access, package installs, or filesystem and network access.
For sandbox setup, providers, and file transfer APIs, see [Sandboxes](/oss/python/deepagents/sandboxes). For the QuickJS runtime and programmatic tool calling, see [Interpreters](/oss/python/deepagents/interpreters).

## [​](#context-management)Context management

The context management component controls what the agent knows, how long it can operate within token limits, and what it retains across sessions. It has four layers:

{' ' * (self.list_depth - 1)}- **[Skills](#skills)**—on-demand domain knowledge loaded progressively from skill files

{' ' * (self.list_depth - 1)}- **[Memory](#memory)**—persistent instructions and preferences loaded at startup from `AGENTS.md` files

{' ' * (self.list_depth - 1)}- **[Summarization and context offloading](#summarization-and-context-offloading)**—automatic compression of conversation history and large tool results

{' ' * (self.list_depth - 1)}- **[Prompt caching](#prompt-caching)**—static prompt sections are cache-eligible to speed up inference and reduce cost on supported models

### [​](#skills)Skills

The harness supports skills that provide specialized workflows and domain knowledge to your deep agent.
**How it works:**

{' ' * (self.list_depth - 1)}- Skills follow the [Agent Skills standard](https://agentskills.io/)

{' ' * (self.list_depth - 1)}- Each skill is a directory containing a `SKILL.md` file with instructions and metadata

{' ' * (self.list_depth - 1)}- Skills can include additional scripts, reference docs, templates, and other resources

{' ' * (self.list_depth - 1)}- Skills use progressive disclosure—they are only loaded when the agent determines they’re useful for the current task

{' ' * (self.list_depth - 1)}- Agent reads frontmatter from each `SKILL.md` file at startup, then reviews full skill content when needed

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Reduces token usage by only loading relevant skills when needed

{' ' * (self.list_depth - 1)}- Bundles capabilities together into larger actions with additional context

{' ' * (self.list_depth - 1)}- Provides specialized expertise without cluttering the system prompt

{' ' * (self.list_depth - 1)}- Enables modular, reusable agent capabilities

For more information, see [Skills](/oss/python/deepagents/skills).

### [​](#memory)Memory

The harness supports persistent memory files that provide extra context to your deep agent across conversations.
These files often contain general coding style, preferences, conventions, and guidelines that help the agent understand how to work with your codebase and follow your preferences.
**How it works:**

{' ' * (self.list_depth - 1)}- Uses [`AGENTS.md` files](https://agents.md/) to provide persistent context

{' ' * (self.list_depth - 1)}- Memory files are always loaded (unlike skills, which use progressive disclosure)

{' ' * (self.list_depth - 1)}- Pass one or more file paths to the `memory` parameter when creating your agent

{' ' * (self.list_depth - 1)}- Files are stored in the agent’s backend (StateBackend, StoreBackend, or FilesystemBackend)

{' ' * (self.list_depth - 1)}- The agent can update memory based on your interactions, feedback, and identified patterns

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Provides persistent context that does not need to be re-specified each conversation

{' ' * (self.list_depth - 1)}- Useful for storing user preferences, project guidelines, or domain knowledge

{' ' * (self.list_depth - 1)}- Always available to the agent, ensuring consistent behavior

For configuration details and examples, see [Memory](/oss/python/deepagents/customization#memory).

### [​](#summarization-and-context-offloading)Summarization and context offloading

The harness manages context so deep agents can handle long-running tasks within token limits while retaining the information they need.
**How it works:**

{' ' * (self.list_depth - 1)}- **Input context**—System prompt, memory, skills, and tool prompts shape what the agent knows at startup

{' ' * (self.list_depth - 1)}- **Compression**—Built-in offloading and summarization keep context within window limits as tasks progress

{' ' * (self.list_depth - 1)}- **Isolation**—Subagents quarantine heavy work and return only results (see [Delegation](#delegation))

{' ' * (self.list_depth - 1)}- **Long-term memory**—Persistent storage across threads via the virtual filesystem

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Enables multi-step tasks that exceed a single context window

{' ' * (self.list_depth - 1)}- Keeps the most relevant information in scope without manual trimming

{' ' * (self.list_depth - 1)}- Reduces token usage through automatic summarization and offloading

For configuration details, see [Context engineering](/oss/python/deepagents/context-engineering).

### [​](#prompt-caching)Prompt caching

For Anthropic models, `create_deep_agent` automatically applies prompt caching to static sections of the system prompt—the base agent instructions, memory, and skill content that repeat on every turn. This avoids reprocessing the same tokens across calls, reducing both latency and cost on long-running agents.
Prompt caching is enabled by default when using an Anthropic model. No configuration is required.
For other providers, see [Middleware integrations](/oss/python/integrations/middleware#official-integrations) for available provider-specific caching middleware.

## [​](#delegation)Delegation

The delegation component enables agents to break large problems into smaller, parallelizable units of work. It has two layers:

{' ' * (self.list_depth - 1)}- **[Task planning](#task-planning)**: a built-in `write_todos` tool for structured task tracking

{' ' * (self.list_depth - 1)}- **[Subagents](#subagents)**: ephemeral child agents that handle isolated subtasks

### [​](#task-planning)Task planning

The harness provides a `write_todos` tool that agents can use to maintain a structured task list.
**Features:**

{' ' * (self.list_depth - 1)}- Track multiple tasks with statuses (`'pending'`, `'in_progress'`, `'completed'`)

{' ' * (self.list_depth - 1)}- Persisted in agent state

{' ' * (self.list_depth - 1)}- Helps agent organize complex multi-step work

{' ' * (self.list_depth - 1)}- Useful for long-running tasks and planning

### [​](#subagents)Subagents

The harness allows the main agent to create ephemeral “subagents” for isolated multi-step tasks.
**Why it’s useful:**

{' ' * (self.list_depth - 1)}- **Context isolation**—Subagent’s work does not clutter main agent’s context

{' ' * (self.list_depth - 1)}- **Parallel execution**—Multiple subagents can run concurrently

{' ' * (self.list_depth - 1)}- **Specialization**—Subagents can have different tools and configurations

{' ' * (self.list_depth - 1)}- **Token efficiency**—Large subtask context is compressed into a single result

**How it works:**

{' ' * (self.list_depth - 1)}- Main agent has a `task` tool

{' ' * (self.list_depth - 1)}- When invoked, it creates a fresh agent instance with its own context

{' ' * (self.list_depth - 1)}- Subagent executes autonomously until completion

{' ' * (self.list_depth - 1)}- Returns a single final report to the main agent

{' ' * (self.list_depth - 1)}- Can use [default `general-purpose` subagent](/oss/python/deepagents/subagents#default-subagent) (enabled by default) or add [custom subagents](/oss/python/deepagents/subagents#custom-subagents)

{' ' * (self.list_depth - 1)}- Subagents are stateless (cannot send multiple messages back)

Running without subagents (no `task` tool)To run an agent without the `task` tool, see [Running without subagents](/oss/python/deepagents/subagents#running-without-subagents). Do not try removing `SubAgentMiddleware` via `excluded_middleware`—that is intentionally rejected. Instead, disable the auto-added subagent via the [harness profile](#harness-profiles) and pass no synchronous subagents via `subagents=`. Async subagents are unaffected.
For more information, see [Subagents](/oss/python/deepagents/subagents).

## [​](#steering)Steering

The steering component gives humans control over agent behavior at runtime.

### [​](#human-in-the-loop)Human-in-the-loop

The harness can pause agent execution at specified tool calls to allow human approval or modification. This feature is opt-in via the `interrupt_on` parameter.
**Configuration:**

{' ' * (self.list_depth - 1)}- Pass `interrupt_on` to `create_deep_agent` with a mapping of tool names to interrupt configurations

{' ' * (self.list_depth - 1)}- Example: `interrupt_on={"edit_file": True}` pauses before every edit

{' ' * (self.list_depth - 1)}- You can provide approval messages or modify tool inputs when prompted

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Safety gates for destructive operations

{' ' * (self.list_depth - 1)}- User verification before expensive API calls

{' ' * (self.list_depth - 1)}- Interactive debugging and guidance

For more information, see [Human-in-the-loop](/oss/python/deepagents/human-in-the-loop).

## [​](#harness-profiles)Harness profiles

The harness can apply a declarative configuration bundle (a `HarnessProfile`) whenever a given provider or model is selected. Profiles tune runtime behavior after the model is built, without requiring per-agent setup code.
**How it works:**

{' ' * (self.list_depth - 1)}- Register a profile under a provider name (`"openai"`) or a `provider:model` key (`"openai:gpt-5.4"`)

{' ' * (self.list_depth - 1)}- `create_deep_agent` looks up and applies the profile when resolving the model

{' ' * (self.list_depth - 1)}- Provider-level and model-level profiles merge at resolution time

**Why it’s useful:**

{' ' * (self.list_depth - 1)}- Package per-provider or per-model defaults (system-prompt tweaks, tool overrides, middleware) in one place

{' ' * (self.list_depth - 1)}- Keep the `create_deep_agent` call site unchanged when switching models

{' ' * (self.list_depth - 1)}- Ship reusable profiles as plugins via entry points

For the full field list, merge semantics, and plugin packaging, see [Profiles](/oss/python/deepagents/profiles).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/harness.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Going to productionPrevious](/oss/python/deepagents/going-to-production)[ModelsNext](/oss/python/deepagents/models)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
