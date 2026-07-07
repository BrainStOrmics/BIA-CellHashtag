# API Documentation Index

Extracted API documentation for LangGraph, LangChain, and Deep Agents.

## Projects

| Project | Summary | Pages |
|---------|---------|-------|
| [LangGraph](langgraph/SUMMARY.md) | Low-level orchestration runtime for building stateful agents | 36 pages |
| [LangChain](langchain/SUMMARY.md) | Agent framework with composable model, tools, prompt, middleware | 71 pages |
| [Deep Agents](deepagent/SUMMARY.md) | Batteries-included agent harness with planning, filesystem, subagents, memory | 57 pages |

## Relationship

```
Deep Agents (agent harness)
  └── built on LangChain (framework)
        └── built on LangGraph (orchestration runtime)
```

- **Deep Agents**: Planning, subagents, filesystem tools, context management on top of LangGraph
- **LangChain**: Abstractions and integrations for models, tools, and agent loops
- **LangGraph**: Durable execution, streaming, human-in-the-loop, and persistence

## Raw HTML Sources

All raw HTML pages are preserved in:
- `raw/langgraph/` — 36 LangGraph pages
- `raw/langchain/` — 71 LangChain pages
- `raw/deepagent/` — 57 DeepAgent pages
