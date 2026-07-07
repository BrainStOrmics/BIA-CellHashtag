# LangSmith Observability

[Deploy with LangSmith](/oss/python/langchain/deploy)

# LangSmith Observability
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.As you build and run agents with LangChain, you need visibility into how they behave: which [tools](/oss/python/langchain/tools) they call, what prompts they generate, and how they make decisions. LangChain agents built with [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) automatically support tracing through [LangSmith](/langsmith/home), a platform for capturing, debugging, evaluating, and monitoring LLM application behavior.
[*Traces*](/langsmith/observability-concepts#traces) record every step of your agent’s execution, from the initial user input to the final response, including all tool calls, model interactions, and decision points. This execution data helps you debug issues, evaluate performance across different inputs, and monitor usage patterns in production.
This guide shows you how to enable tracing for your LangChain agents and use LangSmith to analyze their execution.

## [​](#prerequisites)Prerequisites

Before you begin, ensure you have the following:

{' ' * (self.list_depth - 1)}- **A LangSmith account**: Sign up (for free) or log in at [smith.langchain.com](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-observability).

{' ' * (self.list_depth - 1)}- **A LangSmith API key**: Follow the [Create an API key](/langsmith/create-account-api-key) guide.

## [​](#enable-tracing)Enable tracing

All LangChain agents automatically support LangSmith tracing. To enable it, set the following environment variables:

```
export LANGSMITH_TRACING=true
export LANGSMITH_API_KEY=<your-api-key>

```

## [​](#quickstart)Quickstart

No extra code is needed to log a trace to LangSmith. Just run your agent code as you normally would:

```
from langchain.agents import create_agent

def send_email(to: str, subject: str, body: str):
 """Send an email to a recipient."""
 # ... email sending logic
 return f"Email sent to {to}"

def search_web(query: str):
 """Search the web for information."""
 # ... web search logic
 return f"Search results for: {query}"

agent = create_agent(
 model="gpt-5.4",
 tools=[send_email, search_web],
 system_prompt="You are a helpful assistant that can send emails and search the web."
)

# Run the agent - all steps will be traced automatically
response = agent.invoke({
 "messages": [{"role": "user", "content": "Search for the latest AI news and email a summary to john@example.com"}]
})

```

By default, the trace will be logged to the project with the name `default`. To configure a custom project name, see [Log to a project](/langsmith/log-traces-to-project).

## [​](#trace-selectively)Trace selectively

You may opt to trace specific invocations or parts of your application using LangSmith’s `tracing_context` context manager:

```
import langsmith as ls

# This WILL be traced
with ls.tracing_context(enabled=True):
 agent.invoke({"messages": [{"role": "user", "content": "Send a test email to alice@example.com"}]})

# This will NOT be traced (if LANGSMITH_TRACING is not set)
agent.invoke({"messages": [{"role": "user", "content": "Send another email"}]})

```

## [​](#log-to-a-project)Log to a project

StaticallyYou can set a custom project name for your entire application by setting the `LANGSMITH_PROJECT` environment variable:
```
export LANGSMITH_PROJECT=my-agent-project

```

DynamicallyYou can set the project name programmatically for specific operations:
```
import langsmith as ls

with ls.tracing_context(project_name="email-agent-test", enabled=True):
 response = agent.invoke({
 "messages": [{"role": "user", "content": "Send a welcome email"}]
 })

```

## [​](#add-metadata-to-traces)Add metadata to traces

You can annotate your traces with custom metadata and tags:

```
response = agent.invoke(
 {"messages": [{"role": "user", "content": "Send a welcome email"}]},
 config={
 "tags": ["production", "email-assistant", "v1.0"],
 "metadata": {
 "user_id": "user_123",
 "session_id": "session_456",
 "environment": "production"
 }
 }
)

```

`tracing_context` also accepts tags and metadata for fine-grained control:

```
with ls.tracing_context(
 project_name="email-agent-test",
 enabled=True,
 tags=["production", "email-assistant", "v1.0"],
 metadata={"user_id": "user_123", "session_id": "session_456", "environment": "production"}):
 response = agent.invoke(
 {"messages": [{"role": "user", "content": "Send a welcome email"}]}
 )

```

This custom metadata and tags will be attached to the trace in LangSmith.
To learn more about how to use traces to debug, evaluate, and monitor your agents, see the [LangSmith documentation](/langsmith/home).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/observability.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[LangSmith DeploymentPrevious](/oss/python/langchain/deploy)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
