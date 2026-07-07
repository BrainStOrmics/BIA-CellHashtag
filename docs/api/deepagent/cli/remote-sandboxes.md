# Use remote sandboxes

[Deep Agents Code](/oss/python/deepagents/code/overview)

# Use remote sandboxes
Copy page

Run Deep Agents Code tool execution in LangSmith, Daytona, Modal, Runloop, or AgentCore sandboxes. Install provider extras, set credentials, and use flags and setup scripts.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Deep Agents Code uses the [sandbox as tool](/oss/python/deepagents/sandboxes#sandbox-as-tool-pattern) pattern: the `dcode` process (LLM loop, memory, tool dispatch) runs on your machine, but agent tool calls (`read_file`, `write_file`, `execute`, etc.) target the remote sandbox, not your local filesystem. To get files into the sandbox, use a [setup script](#setup-scripts) or the provider’s file transfer APIs (see [Working with files](/oss/python/deepagents/sandboxes#working-with-files)).
For a deeper look at sandbox architecture, integration patterns, and security best practices, see [Sandboxes](/oss/python/deepagents/sandboxes).
[](#)

Install provider dependency
{' ' * (self.list_depth - 1)}- LangSmith
{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCoreIncluded by default when installing `deepagents-code`. No extra installation needed.
```
uv tool install deepagents-code --with langchain-daytona

```

```
uv tool install deepagents-code --with langchain-modal

```

```
uv tool install deepagents-code --with langchain-runloop

```

```
uv tool install deepagents-code --with langchain-agentcore-codeinterpreter

```
[](#)

Set provider credentials
{' ' * (self.list_depth - 1)}- LangSmith
{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCore
```
export LANGSMITH_API_KEY="your-key"

```

```
export DAYTONA_API_KEY="your-key"

```

```
modal setup

```

```
export RUNLOOP_API_KEY="your-key"

```

```
export AWS_ACCESS_KEY_ID="your-key"
export AWS_SECRET_ACCESS_KEY="your-secret"
export AWS_SESSION_TOKEN="session-token"
export AWS_REGION="us-west-2"

```
[](#)

Run Deep Agents Code with a sandbox
{' ' * (self.list_depth - 1)}- LangSmith
{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCore
```
dcode --sandbox langsmith

```

```
dcode --sandbox daytona

```

```
dcode --sandbox modal

```

```
dcode --sandbox runloop

```

```
dcode --sandbox agentcore

```

## [​](#sandbox-flags-and-examples)Sandbox flags and examples

| 
 | Flag | Description
 | `--sandbox TYPE` | Sandbox provider to use: `langsmith`, `agentcore`, `modal`, `daytona`, or `runloop` (default: `none`)
 | `--sandbox-id ID` | Reuse an existing sandbox by ID instead of creating a new one. Skips creation and cleanup. Refer to your sandbox documentation for more
 | `--sandbox-setup PATH` | Path to a setup script to run inside the sandbox upon creation
Examples:

```
# Create a new Daytona sandbox
dcode --sandbox daytona

# Reuse an existing sandbox (skips creation and cleanup)
dcode --sandbox runloop --sandbox-id dbx_abc123

# Run a setup script after sandbox creation
dcode --sandbox modal --sandbox-setup ./setup.sh

```

## [​](#setup-scripts)Setup scripts

Use `--sandbox-setup` to run a shell script inside the sandbox after creation. This is useful for cloning repos, installing dependencies, and configuring environment variables.
setup.sh
```
#!/bin/bash
set -e

# Clone repository using GitHub token
git clone https://x-access-token:${GITHUB_TOKEN}@github.com/username/repo.git $HOME/workspace
cd $HOME/workspace

# Make environment variables persistent
cat >> ~/.bashrc <<'EOF'
export GITHUB_TOKEN="${GITHUB_TOKEN}"
export OPENAI_API_KEY="${OPENAI_API_KEY}"
cd $HOME/workspace
EOF
source ~/.bashrc

```

Deep Agents Code expands `${VAR}` references in setup scripts using your local environment variables. Store secrets in a local `.env` file for the setup script to access.
Sandboxes isolate code execution, but agents remain vulnerable to prompt injection with untrusted inputs. Use human-in-the-loop approval, short-lived secrets, and trusted setup scripts only. See [Security considerations](/oss/python/deepagents/sandboxes#security-considerations) for details.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/remote-sandboxes.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Memory and SkillsPrevious](/oss/python/deepagents/code/memory-and-skills)[Use subagents in Deep Agents CodeNext](/oss/python/deepagents/code/subagents)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
