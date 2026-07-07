# LangSmith Observability

[Production](/oss/python/langgraph/application-structure)

# LangSmith Observability
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Traces are a series of steps that your application takes to go from input to output. Each of these individual steps is represented by a run. You can use [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-observability) to visualize these execution steps. To use it, [enable tracing for your application](/langsmith/trace-with-langgraph). This enables you to do the following:

{' ' * (self.list_depth - 1)}- [Debug a locally running application](/langsmith/observability-studio#debug-langsmith-traces).

{' ' * (self.list_depth - 1)}- [Evaluate the application performance](/oss/python/langchain/test/evals).

{' ' * (self.list_depth - 1)}- [Monitor the application](/langsmith/dashboards).

## [​](#prerequisites)Prerequisites

Before you begin, ensure you have the following:

{' ' * (self.list_depth - 1)}- **A LangSmith account**: Sign up (for free) or log in at [smith.langchain.com](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-observability).

{' ' * (self.list_depth - 1)}- **A LangSmith API key**: Follow the [Create an API key](/langsmith/create-account-api-key) guide.

## [​](#enable-tracing)Enable tracing

To enable tracing for your application, set the following environment variables:

```
export LANGSMITH_TRACING=true
export LANGSMITH_API_KEY=<your-api-key>

```

By default, the trace will be logged to the project with the name `default`. To configure a custom project name, see [Log to a project](#log-to-a-project).
For more information, see [Trace with LangGraph](/langsmith/trace-with-langgraph).

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

## [​](#use-anonymizers-to-prevent-logging-of-sensitive-data-in-traces)Use anonymizers to prevent logging of sensitive data in traces

You may want to mask sensitive data to prevent it from being logged to LangSmith. You can create [anonymizers](/langsmith/mask-inputs-outputs#rule-based-masking-of-inputs-and-outputs) and apply them to
your graph using configuration. This example will redact anything matching the Social Security Number format XXX-XX-XXXX from traces sent to LangSmith.
Python
```
from langchain_core.tracers.langchain import LangChainTracer
from langgraph.graph import StateGraph, MessagesState
from langsmith import Client
from langsmith.anonymizer import create_anonymizer

anonymizer = create_anonymizer([
 # Matches SSNs
 { "pattern": r"\b\d{3}-?\d{2}-?\d{4}\b", "replace": "<ssn>" }
])

tracer_client = Client(anonymizer=anonymizer)
tracer = LangChainTracer(client=tracer_client)
# Define the graph
graph = (
 StateGraph(MessagesState)
 ...
 .compile()
 .with_config({'callbacks': [tracer]})
)

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/observability.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[LangSmith DeploymentPrevious](/oss/python/langgraph/deploy)[OverviewNext](/oss/python/langgraph/frontend/overview)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
