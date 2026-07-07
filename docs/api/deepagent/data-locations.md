# Data locations

[Deep Agents Code](/oss/python/deepagents/code/overview)

# Data locations
Copy page

Where Deep Agents Code stores configuration, sessions, and customization filesCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Deep Agents Code stores data in two directory hierarchies:

{' ' * (self.list_depth - 1)}- **`~/.deepagents/`** — Deep Agents-specific data (agent memory, skills, sessions)

{' ' * (self.list_depth - 1)}- **`~/.agents/`** — Tool-agnostic data (skills shared across AI CLI tools)

## [​](#directory-structure)Directory structure

```
~/.deepagents/
├── .state/ # Per-machine Deep Agents Code state (managed automatically)
│ ├── sessions.db # SQLite database for conversation checkpoints
│ ├── history.jsonl # Command input history
│ ├── ... # Other markers & credentials
└── {agent}/ # Per-agent directory (default: "agent")
 ├── AGENTS.md # User customizations to agent instructions
 ├── skills/ # User-level skills
 │ └── {skill-name}/
 │ └── SKILL.md
 └── agents/ # Custom subagent definitions
 └── {subagent-name}/
 └── AGENTS.md

~/.agents/ # Tool-agnostic alias (shared across AI CLIs)
└── skills/ # Skills available to any compatible tool
 └── {skill-name}/
 └── SKILL.md

{project}/ # Project-level (in git repo root)
├── AGENTS.md # Project instructions (root-level)
└── .deepagents/
│ ├── AGENTS.md # Project instructions (preferred location)
│ ├── skills/ # Project-specific skills
│ │ └── {skill-name}/
│ │ └── SKILL.md
│ └── agents/ # Project-specific subagents
│ └── {subagent-name}/
│ └── AGENTS.md
└── .agents/ # Tool-agnostic project skills
 └── skills/
 └── {skill-name}/
 └── SKILL.md

```

## [​](#what-goes-where)What goes where

| 
 | Data | Location | Read/Write | Notes
 | **Sessions** | `~/.deepagents/.state/sessions.db` | R/W | SQLite checkpoint database
 | **Input history** | `~/.deepagents/.state/history.jsonl` | R/W | JSON-lines, up/down arrow recall
 | **Base instructions** | Package `default_agent_prompt.md` | R | Immutable, updated with Deep Agents Code upgrades
 | **User customizations** | `~/.deepagents/{agent}/AGENTS.md` | R/W | Appended to base instructions
 | **Project instructions** | `.deepagents/AGENTS.md` or `AGENTS.md` | R | Both loaded if present
 | **User skills** | `~/.deepagents/{agent}/skills/` | R/W | Agent-specific skills
 | **Shared skills** | `~/.agents/skills/` | R | Tool-agnostic, cross-CLI
 | **Project skills** | `.deepagents/skills/` or `.agents/skills/` | R | Project-scoped
 | **Custom subagents** | `~/.deepagents/{agent}/agents/` | R/W | User-defined subagents
 | **Project subagents** | `.deepagents/agents/` | R | Project-defined subagents

## [​](#precedence-rules)Precedence rules

When the same item exists in multiple locations, **higher precedence wins completely** (no merging).

### [​](#skills)Skills

Precedence order (lowest to highest):

{' ' * (self.list_depth - 1)}- `~/.deepagents/{agent}/skills/` — User Deep Agents Code

{' ' * (self.list_depth - 1)}- `~/.agents/skills/` — User tool-agnostic

{' ' * (self.list_depth - 1)}- `.deepagents/skills/` — Project Deep Agents Code

{' ' * (self.list_depth - 1)}- `.agents/skills/` — Project tool-agnostic *(highest)*

When a skill is loaded, Deep Agents Code verifies that the resolved file path stays within one of these directories. Symlinks that resolve outside all skill roots are rejected. To allow symlink targets in additional directories, see [`[skills].extra_allowed_dirs`](/oss/python/deepagents/code/configuration#skills-extra-allowed-directories).

### [​](#subagents)Subagents

Precedence order (lowest to highest):

{' ' * (self.list_depth - 1)}- `~/.deepagents/{agent}/agents/` — User-level

{' ' * (self.list_depth - 1)}- `.deepagents/agents/` — Project-level *(highest)*

Each subagent is an `AGENTS.md` file with YAML frontmatter (`name`, `description`, optional `model`) and a markdown body for the system prompt. See [Use subagents in Deep Agents Code](/oss/python/deepagents/code/subagents) for the full format reference.

### [​](#instructions)Instructions

All instruction sources are **combined** (not overridden):

{' ' * (self.list_depth - 1)}- Package base prompt *(always loaded)*

{' ' * (self.list_depth - 1)}- `~/.deepagents/{agent}/AGENTS.md` *(appended)*

{' ' * (self.list_depth - 1)}- `.deepagents/AGENTS.md` *(appended)*

{' ' * (self.list_depth - 1)}- `AGENTS.md` at project root *(appended)*

## [​](#deepagents-vs-agents)`.deepagents` vs `.agents`

| 
 | Directory | Purpose | When to use
 | `.deepagents/` | Deep Agents Code-specific | Skills and config that use Deep Agents Code-specific features
 | `.agents/` | Tool-agnostic | Skills you want to share across different AI CLI tools
Use `.agents/skills/` for skills that work with any AI coding assistant.
Use `.deepagents/skills/` for skills that rely on Deep Agents-specific tools or conventions.

## [​](#cleaning-up)Cleaning up

| 
 | Need | Action
 | Reset all data | `rm -rf ~/.deepagents`
 | Clear sessions only | `rm ~/.deepagents/.state/sessions.db*`
 | Clear input history | `rm ~/.deepagents/.state/history.jsonl`
 | Clear stored API keys | `rm ~/.deepagents/.state/auth.json`
 | Clear MCP OAuth tokens | `rm -rf ~/.deepagents/.state/mcp-tokens`
 | Re-run first-run onboarding | `rm ~/.deepagents/.state/onboarding_complete`
 | Reset agent instructions | `dcode agents reset --agent {name}`
 | Remove a skill | `rm -rf ~/.deepagents/{agent}/skills/{skill-name}`
Deleting `~/.deepagents/.state/sessions.db` will remove all conversation history and checkpoints.This cannot be undone unless you have a backup of the `sessions.db` file.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/data-locations.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[MCP toolsPrevious](/oss/python/deepagents/code/mcp-tools)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
