# Deep Agents API Summary

Deep Agents is a "batteries-included" agent harness built on LangChain and LangGraph. Provides built-in capabilities for task planning, file system context management, subagent spawning, and long-term memory.

## Key Characteristics

- **Agent harness**: Same core tool-calling loop as other frameworks, but with built-in reliability features
- **Standalone**: `deepagents` library built on LangChain components and LangGraph runtime
- **Context-rich**: Manages large context via filesystem tools and summarization
- **Multi-model**: Supports Google, OpenAI, Anthropic, OpenRouter, Fireworks, Baseten, Ollama

## Core API: `create_deep_agent`

```python
from deepagents import create_deep_agent

agent = create_deep_agent(
    model="anthropic:claude-sonnet-4-6",
    tools=[my_tool],
    system_prompt="You are a helpful assistant",
)
agent.invoke({"messages": [{"role": "user", "content": "..."}]})
```

## Core Capabilities

### Planning
- Built-in `write_todos` tool for task decomposition and progress tracking
- Agents plan work, track status, and revise plans dynamically

### Context Management
- **File system tools**: `ls`, `read_file`, `write_file`, `edit_file`
- **Summarization**: Compacts older conversation history during long runs
- **Data locations**: Store files, intermediate outputs, large results outside context window

### Execution
- **Sandboxes** (`/sandboxes`): Isolated environments with `execute` tool for shell commands
- **Interpreters** (`/interpreters`): QuickJS runtime for tool composition and data transforms
- **MCP tools** (`/mcp`): Model Context Protocol tool integration

### Subagents (`/subagents`)
- Built-in `task` tool for spawning isolated sub-agents
- Async subagents for parallel execution
- Context isolation — subagents don't pollute main agent's context

### Memory & Skills
- **Long-term memory** (`/memory`): Persist information across threads
- **Skills** (`/skills`): Package reusable procedures and domain knowledge

### Control & Safety
- **Backends** (`/backends`): Pluggable storage (in-memory, local disk, durable stores, sandboxes)
- **Permissions** (`/permissions`): Declarative rules restricting filesystem access
- **Human-in-the-loop** (`/human-in-the-loop`): Pause for human approval on sensitive operations

## Usage Modes

### CLI
- **Overview** (`/cli/overview`): Command-line interface for deep agents
- **Configuration** (`/cli/configuration`): CLI config options
- **Subagents** (`/cli/subagents`): CLI subagent management
- **MCP tools** (`/cli/mcp-tools`): CLI tool integration

### Code
- **Overview** (`/code/overview`): Terminal coding agent
- **Configuration** (`/code/configuration`): Code mode config
- **Memory & skills** (`/code/memory-and-skills`): Knowledge persistence
- **Remote sandboxes** (`/code/remote-sandboxes`): Remote execution

### ACP (Agent Client Protocol)
- Integration with code editors like Zed (`/acp`)

## Specialized Use Cases

- **Deep research** (`/deep-research`): Research agent patterns
- **Data analysis** (`/data-analysis`): Data exploration and visualization
- **Content builder** (`/content-builder`): Content generation workflows

## Deployment

- **Deploy overview** (`/deploy/overview`): Production deployment
- **Going to production** (`/going-to-production`): Production readiness guide
- **Event streaming** (`/event-streaming`): Streaming events in production
- **Frontend** (`/frontend/overview`): Frontend integration patterns

## Comparison

- **Deep Agents vs Claude Agent SDK** (`/comparison`): Side-by-side comparison with Anthropic's harness

## Installation

```bash
pip install -U deepagents
```

## Full Documentation Pages

All 57 pages are available in `docs/api/deepagent/` (extracted content) and `docs/api/raw/deepagent/` (raw HTML).
