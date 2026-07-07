# Agent Chat UI

[Production](/oss/python/langgraph/application-structure)

# Agent Chat UI
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[Agent Chat UI](https://github.com/langchain-ai/agent-chat-ui) is a Next.js application that provides a conversational interface for interacting with any LangChain agent. It supports real-time chat, tool visualization, and advanced features like time-travel debugging and state forking. Agent Chat UI works seamlessly with agents created using [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) and provides interactive experiences for your agents with minimal setup, whether you’re running locally or in a deployed context (such as [LangSmith](/langsmith/home)).
Agent Chat UI is open source and can be adapted to your application needs.

You can use generative UI in the Agent Chat UI. For more information, see [Implement generative user interfaces with LangGraph](/langsmith/generative-ui-react).

### [​](#quick-start)Quick start

The fastest way to get started is using the hosted version:

{' ' * (self.list_depth - 1)}- **Visit [Agent Chat UI](https://agentchat.vercel.app)**

{' ' * (self.list_depth - 1)}- **Connect your agent** by entering your deployment URL or local server address

{' ' * (self.list_depth - 1)}- **Start chatting** - the UI will automatically detect and render tool calls and interrupts

### [​](#local-development)Local development

For customization or local development, you can run Agent Chat UI locally:
Use npxClone repository
```
# Create a new Agent Chat UI project
npx create-agent-chat-app --project-name my-chat-ui
cd my-chat-ui

# Install dependencies and start
pnpm install
pnpm dev

```

### [​](#connect-to-your-agent)Connect to your agent

Agent Chat UI can connect to both [local](/oss/python/langgraph/studio#set-up-local-agent-server) and [deployed agents](/oss/python/langgraph/deploy).
After starting Agent Chat UI, you’ll need to configure it to connect to your agent:

{' ' * (self.list_depth - 1)}- **Graph ID**: Enter your graph name (find this under `graphs` in your `langgraph.json` file)

{' ' * (self.list_depth - 1)}- **Deployment URL**: Your Agent server’s endpoint (e.g., `http://localhost:2024` for local development, or your deployed agent’s URL)

{' ' * (self.list_depth - 1)}- **LangSmith API key (optional)**: Add your LangSmith API key (not required if you’re using a local Agent server)

Once configured, Agent Chat UI will automatically fetch and display any interrupted threads from your agent.
Agent Chat UI has out-of-the-box support for rendering tool calls and tool result messages. To customize what messages are shown, see [Hiding Messages in the Chat](https://github.com/langchain-ai/agent-chat-ui?tab=readme-ov-file#hiding-messages-in-the-chat).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/ui.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[LangSmith StudioPrevious](/oss/python/langgraph/studio)[LangSmith DeploymentNext](/oss/python/langgraph/deploy)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
