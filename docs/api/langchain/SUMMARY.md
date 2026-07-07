# LangChain API Summary

LangChain is an agent framework providing `create_agent` — a minimal, highly configurable harness composed of model, tools, prompt, and middleware. Built on top of LangGraph for durable execution, human-in-the-loop, and persistence.

## Key Characteristics

- **Agent = Model + Harness**: Minimal, composable agent architecture
- **Highly configurable**: Extend with middleware for custom behavior
- **Provider-agnostic**: Supports OpenAI, Anthropic, Google, and many more
- **Built on LangGraph**: Inherits durable execution, streaming, persistence

## Core API: `create_agent`

```python
from langchain.agents import create_agent

agent = create_agent(
    model="openai:gpt-5.4",
    tools=[my_tool],
    system_prompt="You are a helpful assistant",
)
result = agent.invoke({"messages": [{"role": "user", "content": "..."}]})
```

## Key Components

### Models (`/models`)
- Standardized interface across providers (OpenAI, Anthropic, Google, etc.)
- Swap providers without changing application code

### Tools (`/tools`)
- Tool definition and composition patterns
- Integration with external APIs and services

### Messages (`/messages`)
- Message types and formats for agent conversations
- Structured output handling

### Middleware
- **Overview** (`/middleware/overview`): Middleware layer for agent behavior
- **Built-in** (`/middleware/built-in`): To-do lists, context compression, etc.
- **Custom** (`/middleware/custom`): Building custom middleware

### Memory
- **Short-term memory** (`/short-term-memory`): Working memory within a session
- **Long-term memory** (`/long-term-memory`): Persistent memory across sessions

## Agent Patterns

### Multi-Agent Systems
- **Handoffs** (`/multi-agent/handoffs`): Agent-to-agent task delegation
- **Router** (`/multi-agent/router`): Route requests to specialized agents
- **Skills** (`/multi-agent/skills`): Agent skill registration and discovery
- **Subagents** (`/multi-agent/subagents`): Spawn isolated sub-agents
- **Custom workflow** (`/multi-agent/custom-workflow`): Build custom multi-agent pipelines
- **Supervisor** (`/supervisor`): Supervisor agent pattern

### RAG & Retrieval
- **RAG** (`/rag`): Retrieval-augmented generation patterns
- **Retrieval** (`/retrieval`): Document retrieval and embedding
- **Knowledge base** (`/knowledge-base`): Knowledge base management

### Streaming & Frontend
- **Streaming** (`/streaming/overview`): Stream agent responses
- **Frontend** (`/frontend/overview`): UI integration patterns
- **Integrations**: AI Elements, Assistant UI, CopilotKit, OpenUI

## Production Features

- **Deploy** (`/deploy`): Deployment patterns
- **Observability** (`/observability`): LangSmith integration
- **Evaluation** (`/evals`): Agent evaluation frameworks
- **Testing** (`/test/`): Unit and integration testing
- **Studio** (`/studio`): Visual workflow editor
- **MCP** (`/mcp`): Model Context Protocol support

## Installation

```bash
pip install -U langchain
```

## Full Documentation Pages

All 71 pages are available in `docs/api/langchain/` (extracted content) and `docs/api/raw/langchain/` (raw HTML).
