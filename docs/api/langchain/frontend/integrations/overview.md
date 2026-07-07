# Overview

[Frontend](/oss/python/langchain/frontend/overview)[Integrations](/oss/python/langchain/frontend/integrations/overview)

# Overview
Copy page

Connect useStream to any React UI component library or generative UI frameworkCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[`useStream`](https://reference.langchain.com/javascript/langchain-react/index/useStream) is UI-agnostic. It returns plain reactive state with messages, tool calls, loading flags, and thread history that you wire to any visual layer you choose. These pages show how different libraries integrate with LangChain frontends, each with a different philosophy for building AI chat and generative UI.

## [​](#integrations)Integrations

## CopilotKit
Full AI chat runtime with structured generative UI support. Add a custom CopilotKit endpoint to your LangGraph deployment, then render dynamic component trees in React.

## AI Elements
Composable shadcn/ui-based components for AI chat. Drop in `Conversation`, `Message`, `Tool`, and `Reasoning` and wire them directly to `stream.messages`.

## assistant-ui
Headless React framework with a full runtime layer. Bridge `useStream` to `AssistantRuntimeProvider` via the `useExternalStoreRuntime` adapter.

## OpenUI
Generative UI library that lets the agent produce complete, interactive dashboards in a declarative component DSL. Purpose-built for data-rich, report-style UIs.

## [​](#choosing-a-library)Choosing a library

Each library fits a slightly different integration model. The choice depends on what kind of UI you’re building:

| 
 | | CopilotKit | AI Elements | assistant-ui | OpenUI
 | **Best for** | Full chat runtime plus structured generative UI | Chat with rich message types | Full-featured chat with minimal setup | Generated dashboards and reports
 | **UI style** | CopilotKit chat shell + custom message renderers | Composable shadcn/ui components | Headless slots + default theme | Prebuilt component library with declarative DSL
 | **Customisation** | Custom backend endpoint, agent context, and renderers | Edit source files directly | Override component slots | Theme via CSS custom properties
 | **Streaming UX** | Runtime-managed chat stream with structured assistant payloads | Component-level progressive render | Built-in thread management | Hoisting — shell appears immediately, data fills in
 | **Tool calls** | Via CopilotKit runtime and custom renderers | `Tool` / `ToolHeader` / `ToolOutput` | Custom via message slots | Inline in the generated UI
 | **Agent format** | Structured assistant responses plus optional Markdown | Any `stream.messages` | Any `stream.messages` | Agent outputs openui-lang text
All four work well with LangChain agents, and the latter three also connect directly to [`useStream`](https://reference.langchain.com/javascript/langchain-react/index/useStream). CopilotKit is especially useful when you want a richer runtime layer and a dedicated endpoint that can sit alongside a [LangGraph](/oss/python/langgraph/overview) deployment.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/integrations/overview.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Generative UIPrevious](/oss/python/langchain/frontend/generative-ui)[CopilotKitNext](/oss/python/langchain/frontend/integrations/copilotkit)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
