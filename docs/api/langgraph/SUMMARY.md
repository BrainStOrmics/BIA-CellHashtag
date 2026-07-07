# LangGraph API Summary

LangGraph is a low-level orchestration framework and runtime for building, managing, and deploying long-running, stateful agents. It focuses entirely on agent **orchestration** — durable execution, streaming, human-in-the-loop, and persistence.

## Key Characteristics

- **Low-level**: Does not abstract prompts or architecture; provides infrastructure for any long-running, stateful workflow
- **Agent orchestration focused**: Built for durable execution, streaming, human-in-the-loop workflows
- **Standalone**: Can be used without LangChain, though commonly used together
- **Built on Pregel/Apache Beam**: Inspired by distributed computation frameworks; public API draws from NetworkX

## Core APIs

### Graph API (`langgraph.graph`)
- `StateGraph`: Define agent workflows as state machines with typed state
- `MessagesState`: Built-in state for message-based agents
- `START` / `END`: Special nodes for graph entry and exit points
- `add_node()`, `add_edge()`, `add_conditional_edges()`: Build graph structure
- `compile()`: Compile graph into executable workflow

### Functional API (`langgraph.func`)
- Decorator-based approach for defining graph nodes as functions
- Lower boilerplate for simple workflows

### Runtime (Pregel)
- `Pregel`: The compiled graph executor
- Handles durable execution, checkpointing, streaming
- `invoke()`, `stream()`, `ainvoke()`, `astream()`: Execution methods

## Capabilities

### Persistence & State Management
- **Persistence** (`/persistence`): Checkpoint state across graph executions
- **Durable execution** (`/durable-execution`): Resume workflows from failure points
- **Fault tolerance** (`/fault-tolerance`): Build agents that persist through failures
- **Time travel** (`/use-time-travel`): Inspect and rewind to any checkpointed state

### Streaming & Events
- **Streaming** (`/streaming`): Stream output tokens and intermediate states
- **Event streaming** (`/event-streaming`): Fine-grained event emission during execution

### Human-in-the-Loop
- **Interrupts** (`/interrupts`): Pause execution for human approval
- **Human-in-the-loop** (`/human-in-the-loop`): Full HITL workflow patterns

### Memory & Subgraphs
- **Memory** (`/memory`): Short-term and long-term memory patterns
- **Subgraphs** (`/use-subgraphs`): Nest graphs for modular workflow composition

## Production Features

- **Application structure** (`/application-structure`): Best practices for organizing production code
- **Test** (`/test`): Testing patterns for LangGraph workflows
- **Deploy** (`/deploy`): LangSmith deployment integration
- **Observability** (`/observability`): LangSmith tracing and monitoring
- **Studio** (`/studio`): LangSmith visual workflow editor
- **UI** (`/ui`): Agent chat UI
- **Local server** (`/local-server`): Run LangGraph server locally

## Installation

```bash
pip install -U langgraph
```

## Related Projects

- **LangChain**: Agent framework built on top of LangGraph with higher-level abstractions
- **Deep Agents**: Batteries-included agent harness with planning, filesystem, subagents, memory
- **LangSmith**: Platform for tracing, evaluation, and deployment

## Full Documentation Pages

All 36 pages are available in `docs/api/langgraph/` (extracted content) and `docs/api/raw/langgraph/` (raw HTML).
