# Overview

[Middleware](/oss/python/langchain/middleware/overview)

# Overview
Copy page

Control and customize agent execution at every stepCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Middleware provides a way to more tightly control what happens inside the agent. Middleware is useful for the following:

{' ' * (self.list_depth - 1)}- Tracking agent behavior with logging, analytics, and debugging.

{' ' * (self.list_depth - 1)}- Transforming prompts, [tool selection](/oss/python/langchain/middleware/built-in#llm-tool-selector), and output formatting.

{' ' * (self.list_depth - 1)}- Adding [retries](/oss/python/langchain/middleware/built-in#tool-retry), [fallbacks](/oss/python/langchain/middleware/built-in#model-fallback), and early termination logic.

{' ' * (self.list_depth - 1)}- Applying [rate limits](/oss/python/langchain/middleware/built-in#model-call-limit), guardrails, and [PII detection](/oss/python/langchain/middleware/built-in#pii-detection).

Add middleware by passing them to [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent):

```
from langchain.agents import create_agent
from langchain.agents.middleware import SummarizationMiddleware, HumanInTheLoopMiddleware

agent = create_agent(
 model="gpt-5.4",
 tools=[...],
 middleware=[
 SummarizationMiddleware(...),
 HumanInTheLoopMiddleware(...)
 ],
)

```

## [​](#the-agent-loop)The agent loop

The core agent loop involves calling a model, letting it choose tools to execute, and then finishing when it calls no more tools:

![Core agent loop diagram](https://mintcdn.com/langchain-5e9cc07a/Tazq8zGc0yYUYrDl/oss/images/core_agent_loop.png?fit=max&auto=format&n=Tazq8zGc0yYUYrDl&q=85&s=ac72e48317a9ced68fd1be64e89ec063)
Middleware exposes hooks before and after each of those steps:

![Middleware flow diagram](https://mintcdn.com/langchain-5e9cc07a/RAP6mjwE5G00xYsA/oss/images/middleware_final.png?fit=max&auto=format&n=RAP6mjwE5G00xYsA&q=85&s=eb4404b137edec6f6f0c8ccb8323eaf1)

## [​](#use-middleware-inside-a-langgraph-workflow)Use middleware inside a LangGraph workflow

Middleware is not a separate runtime: hooks run inside the compiled [LangGraph](/oss/python/langgraph/overview) that [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) returns. You can drop the whole agent (middleware and all) into a larger [StateGraph](https://reference.langchain.com/python/langgraph/graph/state/StateGraph) as a node or subgraph, and every middleware hook continues to run.
Reach for this pattern when the surrounding topology is more than a standard “loop until done”: classifying input before routing to one of several agents, fanning out work in parallel, or stitching agent calls together with deterministic steps.
`HumanInTheLoopMiddleware` matches against each tool’s `.name`. In Python, `@tool`-decorated functions take their name from the function (so the key below is `"send_email"`); in TypeScript, the key matches the `name` you pass to `tool({...}, { name })`.

```
from langchain.agents import AgentState, create_agent
from langchain.agents.middleware import HumanInTheLoopMiddleware
from langgraph.graph import START, StateGraph

# Assumes read_email, send_email, classify_node, and route are defined elsewhere.
email_agent = create_agent(
 model="claude-sonnet-4-6",
 tools=[read_email, send_email],
 middleware=[HumanInTheLoopMiddleware(interrupt_on={"send_email": True})],
)

graph = (
 StateGraph(AgentState)
 .add_node("classify", classify_node)
 .add_node("email_agent", email_agent)
 .add_edge(START, "classify")
 .add_conditional_edges("classify", route)
 .compile()
)

```

The HITL interrupt, summarization, PII redaction, retries, and any custom hooks all travel with the agent node. See [Use subgraphs](/oss/python/langgraph/use-subgraphs) for the full set of composition patterns, including subgraph checkpointer scoping (per-invocation versus per-thread).

## [​](#additional-resources)Additional resources

## Built-in middleware
Explore built-in middleware for common use cases.

## Custom middleware
Build your own middleware with hooks and decorators.

## Middleware API reference
Complete API reference for middleware.

## Middleware integrations
Provider-specific middleware for Anthropic, AWS, OpenAI, and more.

## Testing agents
Test your agents with LangSmith.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/middleware/overview.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Structured outputPrevious](/oss/python/langchain/structured-output)[Prebuilt middlewareNext](/oss/python/langchain/middleware/built-in)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
