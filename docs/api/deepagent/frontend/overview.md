# Overview

[Frontend](/oss/python/deepagents/frontend/overview)

# Overview
Copy page

Build UIs that display real-time subagent streams, task progress, and sandbox for Deep AgentsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Build frontends that visualize deep agent workflows in real time. These patterns show how to render subagent progress, task planning, streaming content, and IDE-like sandbox experiences from agents created with `createDeepAgent`.

## [​](#architecture)Architecture

Deep Agents use a coordinator-worker architecture. The main agent plans tasks and delegates to specialized subagents, each running in isolation. On the frontend, `useStream` surfaces both the coordinator’s messages and each subagent’s streaming state.

```
from deepagents import create_deep_agent

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
 subagents=[
 {
 "name": "researcher",
 "description": "Research assistant",
 }
 ],
)

```

On the frontend, connect with `useStream` the same way as with `createAgent`. Deep agent patterns use additional `useStream` features like `stream.subagents`, `stream.values.todos`, and `filterSubagentMessages` to render subagent-specific UIs.

```
import { useStream } from "@langchain/react";

function App() {
 const stream = useStream<typeof agent>({
 apiUrl: "http://localhost:2024",
 assistantId: "agent",
 });

 // Deep agent state beyond messages
 const todos = stream.values?.todos;
 const subagents = stream.subagents;
}

```

## [​](#patterns)Patterns

## Subagent streaming
Display specialist subagents with streaming content, progress tracking, and collapsible cards.

## Todo list
Track agent progress with a real-time todo list synced from agent state.

## Sandbox
Build an IDE-like UI with a file browser, code viewer, and diff panel backed by a sandbox.

## [​](#related-patterns)Related patterns

The [LangChain frontend patterns](/oss/python/langchain/frontend/overview), including
markdown messages, tool calling, and human-in-the-loop, all work with deep
agents too. Deep Agents are built on the same LangGraph runtime, so
`useStream` provides the same core API.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/frontend/overview.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[StreamingPrevious](/oss/python/deepagents/streaming)[Subagent streamingNext](/oss/python/deepagents/frontend/subagent-streaming)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
