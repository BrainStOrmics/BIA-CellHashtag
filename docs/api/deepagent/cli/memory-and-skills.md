# Memory and Skills

[Deep Agents Code](/oss/python/deepagents/code/overview)

# Memory and Skills
Copy page

Persistent memory, AGENTS.md files, and reusable skills for Deep Agents Code, including creation, discovery, and invocation.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.There are two primary ways to customize any agent:

{' ' * (self.list_depth - 1)}- 
**[Memory](#memory)**: `AGENTS.md` files and auto-saved memories that persist across sessions. Use memory for general coding style, preferences, and learned conventions.

{' ' * (self.list_depth - 1)}- 
**[Skills](#skills)**: Global and project-specific context, conventions, guidelines, or instructions. Use skills for context that applies when performing specific tasks.

Use `/remember` to explicitly prompt the agent to update its memory and skills from the current conversation.
Building a custom agent with the SDK? See [Memory](/oss/python/deepagents/memory) for programmatic memory backends.

## [​](#memory)Memory

### [​](#automatic-memory)Automatic memory

As you use the agent, it automatically stores information in `~/.deepagents/<agent_name>/memories/` as markdown files using a memory-first protocol:

{' ' * (self.list_depth - 1)}- **Research**: Searches memory for relevant context before starting tasks

{' ' * (self.list_depth - 1)}- **Response**: Checks memory when uncertain during execution

{' ' * (self.list_depth - 1)}- **Learning**: Automatically saves new information for future sessions

The agent organizes its memories by topic with descriptive filenames:

```
~/.deepagents/backend-dev/memories/
├── api-conventions.md
├── database-schema.md
└── deployment-process.md

```

When you teach the agent conventions:

```
dcode --agent backend-dev
> Our API uses snake_case and includes created_at/updated_at timestamps

```

It remembers for future sessions:

```
> Create a /users endpoint
# Applies conventions without prompting

```

### [​](#agents-md-files)AGENTS.md files

[`AGENTS.md` files](https://agents.md/) provide persistent context that is always loaded at session start:

{' ' * (self.list_depth - 1)}- **Global**: `~/.deepagents/<agent_name>/AGENTS.md` — loaded every session.

{' ' * (self.list_depth - 1)}- **Project**: `.deepagents/AGENTS.md` in any git project root — loaded when Deep Agents Code is run from within that project.

Both files are appended to the system prompt at startup.

### [​](#how-memory-works)How memory works

The agent may also read its memory files when answering project-specific questions or when you reference past work or patterns.
The agent will update `AGENTS.md` as you provide information on how it should behave, feedback on its work, or instructions to remember something.
It will also update its memory if it identifies patterns or preferences from your interactions.
To add more structured project knowledge in additional memory files, add them in `.deepagents/` and reference them in the `AGENTS.md` file.
You must reference additional files in the `AGENTS.md` file for the agent to be aware of them.
The additional files will not be read on startup but the agent can reference and update them when needed.

### [​](#when-to-use-global-vs-project-agents-md)When to use global vs. project AGENTS.md

Use a global `AGENTS.md` (`~/.deepagents/agent/AGENTS.md`) for:

{' ' * (self.list_depth - 1)}- Your personality, style, and universal coding preferences

{' ' * (self.list_depth - 1)}- General tone and communication style

{' ' * (self.list_depth - 1)}- Universal coding preferences (formatting, type hints, etc.)

{' ' * (self.list_depth - 1)}- Tool usage patterns that apply everywhere

{' ' * (self.list_depth - 1)}- Workflows and methodologies that don’t change per-project

Use a project `AGENTS.md` (`.deepagents/AGENTS.md` in project root) for:

{' ' * (self.list_depth - 1)}- Project-specific context and conventions

{' ' * (self.list_depth - 1)}- Project architecture and design patterns

{' ' * (self.list_depth - 1)}- Coding conventions specific to this codebase

{' ' * (self.list_depth - 1)}- Testing strategies and deployment processes

{' ' * (self.list_depth - 1)}- Team guidelines and project structure

## [​](#skills)Skills

Skills are reusable agent capabilities that provide specialized workflows and domain knowledge.
You can use [skills](/oss/python/deepagents/skills) to provide your deep agent with new capabilities and expertise.
Deep agent skills follow the [Agent Skills standard](https://agentskills.io/).
Once you have added skills your deep agent will automatically make use of them and update them as you use the agent and provide it with additional information.

### [​](#add-skills)Add skills

1[](#)

Create a skill
```
# User skill (stored in ~/.deepagents/<agent_name>/skills/)
dcode skills create test-skill

# Project skill (stored in .deepagents/skills/)
dcode skills create test-skill --project

```
This generates:
```
skills/
└── test-skill
 └── SKILL.md

```
2[](#)

Edit SKILL.mdOpen the generated `SKILL.md` and edit the file to include your instructions.3[](#)

Add optional resourcesOptionally add additional scripts or other resources to the `test-skill` folder. For more information, see [Examples](/oss/python/deepagents/skills#example).
You can also copy existing skills directly to the agent’s folder:

```
mkdir -p ~/.deepagents/<agent_name>/skills
cp -r examples/skills/web-research ~/.deepagents/<agent_name>/skills/

```

### [​](#install-community-skills)Install community skills

You can use tools like Vercel’s [Skills CLI](https://github.com/vercel-labs/skills) to install community [Agent Skills](https://agentskills.io/) in your environment and make them available to your deep agents:

```
# Install a skill globally
npx skills add vercel-labs/agent-skills --skill web-design-guidelines -a deepagents -g -y

# List installed skills
npx skills ls -a deepagents -g

```

Global installs (`-g`) symlink skills into `~/.deepagents/agent/skills/` — the default agent’s user-level skills directory. Project-level installs (omit `-g`) place skills in `.deepagents/skills/` relative to the current directory, making them available to any agent running in that project regardless of agent name.
Global installs target the default `agent` directory only. If you use a custom-named agent, either use project-level installs or manually symlink the skill into `~/.deepagents/{your-agent}/skills/`.

### [​](#skill-discovery)Skill discovery

At startup, Deep Agents Code discovers skills from both Deep Agents and shared alias directories:

```
~/.deepagents/<agent_name>/skills/
~/.agents/skills/
.deepagents/skills/
.agents/skills/
~/.claude/skills/ (experimental)
.claude/skills/ (experimental)

```

When duplicate skill names exist, later-precedence directories override earlier ones (see [App data](/oss/python/deepagents/code/data-locations#skills)).
For project-specific skills, the project’s root folder must have a `.git` folder.
When you start Deep Agents Code from anywhere within the project’s folder, it will find the project’s root folder by checking for a containing `.git` folder.
For each skill, Deep Agents Code reads the name and the description from the `SKILL.md` file’s frontmatter.
As you use Deep Agents Code, if a task matches the skill’s description, the agent will read the skill file and follow its instructions.
You can also invoke a skill directly with `/skill:<name> [args]`. Skill discovery runs at startup and again on `/reload`.

### [​](#invoke-a-skill-from-the-command-line)Invoke a skill from the command line

Use `--skill` to invoke a skill at launch without typing a slash command interactively:

```
# Open the TUI and immediately run a skill
dcode --skill code-review

# Pass a request to the skill with -m
dcode --skill code-review -m 'review the auth module'

# Pipe content into a skill
cat diff.txt | dcode --skill code-review

# Pipe content and add a request
cat diff.txt | dcode --skill code-review -m 'focus on security'

```

`--skill` also works in non-interactive mode:

```
# Run a skill headlessly
dcode --skill code-review -n 'review this patch'

# Quiet mode (only agent output on stdout)
dcode --skill code-review -n 'review this patch' -q

```

`--skill` with `--quiet` or `--no-stream` requires `-n` (non-interactive mode).

### [​](#list-skills)List skills

```
# List all user skills
dcode skills list

# List project skills
dcode skills list --project

# Get detailed info about a specific skill
dcode skills info test-skill
dcode skills info test-skill --project

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/memory-and-skills.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Deep Agents CodePrevious](/oss/python/deepagents/code/overview)[Use remote sandboxesNext](/oss/python/deepagents/code/remote-sandboxes)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
