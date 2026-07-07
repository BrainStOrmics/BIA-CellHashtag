# LangChain Documentation

## Summary
See [SUMMARY.md](SUMMARY.md) for comprehensive API overview.

## Get Started
- [Overview](overview.md) — What is LangChain, create_agent harness
- [Install](install.md) — Installation guide
- [Quickstart](quickstart.md) — Build your first agent
- [Academy](academy.md) — LangChain Academy courses
- [Get help](get-help.md) — Support and community
- [Philosophy](philosophy.md) — LangChain design philosophy
- [Component architecture](component-architecture.md) — Architecture overview

## Core Components
- [Models](models.md) — Model providers and standardized interface
- [Tools](tools.md) — Tool definition and composition
- [Messages](messages.md) — Message types and formats
- [Agents](agents/) — create_agent harness and agent loop
- [Runtime](runtime.md) — Agent runtime
- [MCP](mcp.md) — Model Context Protocol

## Middleware
- [Overview](middleware/overview.md) — Middleware layer
- [Built-in](middleware/built-in.md) — Built-in middleware (to-do lists, etc.)
- [Custom](middleware/custom.md) — Custom middleware

## Memory
- [Short-term memory](short-term-memory.md) — Session working memory
- [Long-term memory](long-term-memory.md) — Persistent cross-session memory

## Multi-Agent Patterns
- [Index](multi-agent/index.md) — Multi-agent overview
- [Handoffs](multi-agent/handoffs.md) — Agent-to-agent delegation
- [Handoffs: customer support](multi-agent/handoffs-customer-support.md) — CS example
- [Router](multi-agent/router.md) — Request routing
- [Router: knowledge base](multi-agent/router-knowledge-base.md) — KB routing example
- [Skills](multi-agent/skills.md) — Skill registration/discovery
- [Skills: SQL assistant](multi-agent/skills-sql-assistant.md) — SQL example
- [Subagents](multi-agent/subagents.md) — Spawn isolated sub-agents
- [Subagents: personal assistant](multi-agent/subagents-personal-assistant.md) — PA example
- [Custom workflow](multi-agent/custom-workflow.md) — Custom multi-agent
- [Supervisor](supervisor.md) — Supervisor pattern

## RAG & Retrieval
- [RAG](rag.md) — Retrieval-augmented generation
- [Retrieval](retrieval.md) — Document retrieval
- [Knowledge base](knowledge-base.md) — Knowledge base management

## Streaming
- [Overview](streaming/overview.md) — Streaming overview
- [Frontend](streaming/frontend.md) — Streaming in frontend

## Frontend
- [Overview](frontend/overview.md) — Frontend integration
- [Branching chat](frontend/branching-chat.md) — Branching conversations
- [Generative UI](frontend/generative-ui.md) — Dynamic UI generation
- [Human-in-the-loop](frontend/human-in-the-loop.md) — HITL in frontend
- [Join/rejoin](frontend/join-rejoin.md) — Conversation join patterns
- [Markdown messages](frontend/markdown-messages.md) — Markdown rendering
- [Message queues](frontend/message-queues.md) — Message queue patterns
- [Reasoning tokens](frontend/reasoning-tokens.md) — Reasoning token display
- [Structured output](frontend/structured-output.md) — Structured output in UI
- [Time travel](frontend/time-travel.md) — State rewinding in UI
- [Tool calling](frontend/tool-calling.md) — Tool call visualization
- [Integrations](frontend/integrations.md) — Frontend library integrations
  - [Overview](frontend/integrations/overview.md)
  - [AI Elements](frontend/integrations/ai-elements.md)
  - [Assistant UI](frontend/integrations/assistant-ui.md)
  - [CopilotKit](frontend/integrations/copilotkit.md)
  - [OpenUI](frontend/integrations/openui.md)

## Production
- [Deploy](deploy.md) — Deployment patterns
- [Observability](observability.md) — LangSmith integration
- [Studio](studio.md) — Visual workflow editor
- [UI](ui.md) — Agent chat UI
- [Human-in-the-loop](human-in-the-loop.md) — HITL patterns
- [Harness](harness.md) — Agent harness details
- [Guardrails](guardrails.md) — Safety and guardrails
- [Event streaming](event-streaming.md) — Event emission
- [Structured output](structured-output.md) — Structured response handling
- [Context engineering](context-engineering.md) — Context management
- [Deep agent from scratch](deep-agent-from-scratch.md) — Build deep agent manually
- [Voice agent](voice-agent.md) — Voice interaction patterns
- [SQL agent](sql-agent.md) — SQL agent patterns

## Testing & Evaluation
- [Index](test/index.md) — Testing overview
- [Unit testing](test/unit-testing.md)
- [Integration testing](test/integration-testing.md)
- [Evals](test/evals.md) — Evaluation frameworks
- [Evals](evals.md) — Agent evaluation

## Reference
- [Changelog](changelog-py.md) — Python changelog
