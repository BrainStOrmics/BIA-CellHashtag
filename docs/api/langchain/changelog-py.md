# Changelog

[Releases](/oss/python/versioning)

# Changelog
Copy page

Log of updates and improvements to our Python packagesCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.**Subscribe**: Our changelog includes an [RSS feed](https://docs.langchain.com/oss/python/releases/changelog/rss.xml) that can integrate with [Slack](https://slack.com/help/articles/218688467-Add-RSS-feeds-to-Slack), [email](https://zapier.com/apps/email/integrations/rss/1441/send-new-rss-feed-entries-via-email), Discord bots like [Readybot](https://readybot.io/) or [RSS Feeds to Discord Bot](https://rss.app/en/bots/rssfeeds-discord-bot), and other subscription tools.
[​](#may-12-2026)May 12, 2026deepagents

## [​](#deepagents-v0-6-0)`deepagents` v0.6.0

{' ' * (self.list_depth - 1)}- **[`CodeInterpreterMiddleware`](/oss/python/deepagents/interpreters)**: (experimental) `deepagents` now supports code execution and programmatic tool calling through a scoped QuickJS runtime.

{' ' * (self.list_depth - 1)}- Supports `version="v3"` in `stream_events` / `astream_events`. Refer to the [event streaming](/oss/python/deepagents/event-streaming) guide for details.

[​](#may-12-2026-2)May 12, 2026langchain

## [​](#langchain-v1-3-0)`langchain` v1.3.0
This release adds support for `version="v3"` in `stream_events` / `astream_events` for `langchain` agents. Refer to the [event streaming](/oss/python/langchain/event-streaming) guide for details.
[​](#may-12-2026-3)May 12, 2026langgraph

## [​](#langgraph-v1-2-0)`langgraph` v1.2.0
This release adds finer-grained control over node execution (timeouts, error recovery, and graceful shutdown), a new channel type that cuts checkpoint overhead for long-running threads, and a new content-block-centric streaming API (v3) with typed, per-channel projections.

{' ' * (self.list_depth - 1)}- 
**[`DeltaChannel`](/oss/python/langgraph/pregel#deltachannel-beta) (beta)**: A new channel type that stores only the incremental delta at each step rather than re-serializing the full accumulated value. Most useful for channels that grow large over time, for example a message list in a long-running thread. Use `snapshot_frequency=K` to write a full snapshot every K steps and bound read latency.

{' ' * (self.list_depth - 1)}- 
**[Per-node timeouts](/oss/python/langgraph/fault-tolerance#timeouts)**: Pass `timeout=` to [`add_node`](https://reference.langchain.com/python/langgraph/graph/state/StateGraph/add_node) to cap how long a single attempt may run. Set a hard wall-clock limit (`run_timeout`), an idle limit that resets on progress (`idle_timeout`), or both via [`TimeoutPolicy`](https://reference.langchain.com/python/langgraph/types/TimeoutPolicy). When the limit fires, LangGraph raises [`NodeTimeoutError`](https://reference.langchain.com/python/langgraph/errors/NodeTimeoutError), clears writes from that attempt, and hands off to the retry policy. Async nodes only.

{' ' * (self.list_depth - 1)}- 
**[Node-level error handlers](/oss/python/langgraph/fault-tolerance#error-handling)**: Pass `error_handler=` to [`add_node`](https://reference.langchain.com/python/langgraph/graph/state/StateGraph/add_node) to run a recovery function after all retries are exhausted. The handler receives a typed [`NodeError`](https://reference.langchain.com/python/langgraph/errors/NodeError) and can return a [`Command`](https://reference.langchain.com/python/langgraph/types/Command) to update state and route to a different node, useful for Saga/compensation patterns.

{' ' * (self.list_depth - 1)}- 
**[Graceful shutdown](/oss/python/langgraph/durable-execution#graceful-shutdown)**: Stop an in-flight run cooperatively after the current superstep completes, and save a resumable checkpoint. Create a [`RunControl`](https://reference.langchain.com/python/langgraph/runtime/RunControl) and call `request_drain()` from any thread; the run raises `GraphDrained` and can be resumed later with the same config.

{' ' * (self.list_depth - 1)}- 
**New event streaming API (beta)**: Pass `version="v3"` to `stream_events()` / `astream_events()` for a content-block-centric protocol with typed, per-channel projections (`run.values`, `run.messages`, `run.lifecycle`, `run.subgraphs`) plus opt-in transformers for updates, custom events, checkpoints, tasks, and debug. `run.messages` yields one `ChatModelStream` per LLM call with typed sub-projections for text, reasoning, tool calls, and usage. `version="v1"` and `version="v2"` are unchanged.

Timeouts and error handlers are Python-only; retry policies continue to work in both Python and TypeScript.
[​](#apr-7-2026)Apr 7, 2026deepagents

## [​](#deepagents-v0-5-0)`deepagents` v0.5.0

{' ' * (self.list_depth - 1)}- 
**[Async subagents](/oss/python/deepagents/async-subagents)**: Deep Agents can launch non-blocking background tasks, so users can continue interacting with the agent while subagents work concurrently. Requires [LangSmith Deployment](/langsmith/deployment) for sub-agents.

{' ' * (self.list_depth - 1)}- 
**Multi-modal support**: The `read_file` tool now supports PDFs, audio, and video files in addition to images.

{' ' * (self.list_depth - 1)}- 
**Backend changes**: We’ve made backward-compatible changes to the Deep Agents [backend protocol](https://github.com/langchain-ai/deepagents/blob/main/libs/deepagents/deepagents/backends/protocol.py):

{' ' * (self.list_depth - 1)}- Updated the file format stored in [State and Store backends](/oss/python/deepagents/backends) to support binary files.

{' ' * (self.list_depth - 1)}- Improved error propagation from backends to tools.

{' ' * (self.list_depth - 1)}- You can now instantiate `StateBackend()` and `StoreBackend()` directly. Specifying with a factory (e.g., `backend=(lambda rt: StateBackend(rt))`) is deprecated.

{' ' * (self.list_depth - 1)}- 
**Anthropic prompt caching improvements**: We’ve made some improvements to improve prompt caching performance for Anthropic models.

[​](#mar-10-2026)Mar 10, 2026langgraph

## [​](#langgraph-v1-1-0)`langgraph` v1.1.0

{' ' * (self.list_depth - 1)}- 
**Type-safe streaming (`version="v2"`)**: Pass `version="v2"` to `stream()` / `astream()` for unified `StreamPart` output with `type`, `ns`, and `data` keys on every chunk. Each mode has its own `TypedDict`, all importable from `langgraph.types`. See [streaming docs](/oss/python/langgraph/streaming#stream-output-format-v2).

{' ' * (self.list_depth - 1)}- 
**Type-safe invoke (`version="v2"`)**: Pass `version="v2"` to `invoke()` / `ainvoke()` to get a `GraphOutput` object with `.value` and `.interrupts` attributes. See [invoke docs](/oss/python/langgraph/streaming#v2-invoke-format).

{' ' * (self.list_depth - 1)}- 
**Pydantic and dataclass coercion**: With `version="v2"`, `invoke()` and `values`-mode stream output are automatically coerced to your declared Pydantic model or dataclass type.

{' ' * (self.list_depth - 1)}- 
**Fixed time travel with interrupts and subgraphs**: Replays no longer reuse stale `RESUME` values, and subgraphs correctly restore the checkpoint for the parent’s historical state.

{' ' * (self.list_depth - 1)}- 
**Fully backwards compatible**: `version="v2"` is opt-in. `GraphOutput` supports deprecated dict-style access for gradual migration.

[​](#feb-10-2026)Feb 10, 2026deepagents

## [​](#deepagents-v0-4-0)`deepagents` v0.4.0

{' ' * (self.list_depth - 1)}- New integration packages for pluggable sandboxes: [`langchain-modal`](https://pypi.org/project/langchain-modal/), [`langchain-daytona`](https://pypi.org/project/langchain-daytona/), and [`langchain-runloop`](https://pypi.org/project/langchain-runloop/). See [sandboxes guide](/oss/python/deepagents/sandboxes) and example [data analysis tutorial](/oss/python/deepagents/data-analysis).

{' ' * (self.list_depth - 1)}- Changes to [conversation history summarization](/oss/python/deepagents/context-engineering#summarization):

{' ' * (self.list_depth - 1)}- Summarization now happens in the model node via `wrap_model_call` events. Due to this we retain the full message history in the graph state.

{' ' * (self.list_depth - 1)}- More accurate token counting.

{' ' * (self.list_depth - 1)}- Summarization will now automatically trigger if a chat model raises a [`ContextOverflowError`](https://reference.langchain.com/python/langchain-core/exceptions/ContextOverflowError) (defined in `langchain-core`). Currently `langchain-anthropic` and `langchain-openai` support this.

{' ' * (self.list_depth - 1)}- We now default to the Responses API for model strings prefixed with `"openai:"`.

Disable data retention with the Responses API
```
from langchain.chat_models import init_chat_model

agent = create_deep_agent(
 model=init_chat_model(
 "openai:...",
 use_responses_api=True,
 store=False,
 include=["reasoning.encrypted_content"],
 )
)

```

[​](#dec-15-2025)Dec 15, 2025langchainintegrations

## [​](#langchain-v1-2-0)`langchain` v1.2.0

{' ' * (self.list_depth - 1)}- [`create_agent`](/oss/python/langchain/agents): Simplified support for provider-specific tool parameters and definitions via a new [`extras`](https://reference.langchain.com/python/langchain/tools/#langchain.tools.BaseTool.extras) attribute on [tools](/oss/python/langchain/tools). Examples:

{' ' * (self.list_depth - 1)}- Provider-specific configuration such as Anthropic’s [programmatic tool calling](/oss/python/integrations/chat/anthropic#programmatic-tool-calling) and [tool search](/oss/python/integrations/chat/anthropic#tool-search).

{' ' * (self.list_depth - 1)}- Built-in tools that are executed client-side, as supported by [Anthropic](/oss/python/integrations/chat/anthropic#built-in-tools), [OpenAI](/oss/python/integrations/chat/openai#responses-api), and other providers.

{' ' * (self.list_depth - 1)}- Support for strict schema-adherence in agent `response_format` (see [`ProviderStrategy`](/oss/python/langchain/structured-output#provider-strategy) docs).

[​](#dec-8-2025)Dec 8, 2025langchainintegrations

## [​](#langchain-google-genai-v4-0-0)`langchain-google-genai` v4.0.0
We’ve re-written the Google GenAI integration to use Google’s consolidated Generative AI SDK, which provides access to the Gemini API and Vertex AI Platform under the same interface. This includes minimal breaking changes as well as deprecated packages in `langchain-google-vertexai`.See the full [release notes and migration guide](https://github.com/langchain-ai/langchain-google/discussions/1422) for details.
[​](#nov-25-2025)Nov 25, 2025langchain

## [​](#langchain-v1-1-0)`langchain` v1.1.0

{' ' * (self.list_depth - 1)}- [Model profiles](/oss/python/langchain/models#model-profiles): Chat models now expose supported features and capabilities through a `.profile` attribute. These data are derived from [models.dev](https://models.dev), an open source project providing model capability data.

{' ' * (self.list_depth - 1)}- [Summarization middleware](/oss/python/langchain/middleware/built-in#summarization): Updated to support flexible trigger points using model profiles for context-aware summarization.

{' ' * (self.list_depth - 1)}- [Structured output](/oss/python/langchain/structured-output): `ProviderStrategy` support (native structured output) can now be inferred from model profiles.

{' ' * (self.list_depth - 1)}- [`SystemMessage` for `create_agent`](/oss/python/langchain/middleware/custom#dynamic-prompt): Support for passing `SystemMessage` instances directly to `create_agent`’s `system_prompt` parameter, enabling advanced features like cache control and structured content blocks.

{' ' * (self.list_depth - 1)}- [Model retry middleware](/oss/python/langchain/middleware/built-in#model-retry): New middleware for automatically retrying failed model calls with configurable exponential backoff.

{' ' * (self.list_depth - 1)}- [Content moderation middleware](/oss/python/integrations/middleware/openai#content-moderation): OpenAI content moderation middleware for detecting and handling unsafe content in agent interactions. Supports checking user input, model output, and tool results.

[​](#oct-20-2025)Oct 20, 2025langchainlanggraph

## [​](#v1-0-0)v1.0.0

### [​](#langchain)`langchain`

{' ' * (self.list_depth - 1)}- [Release notes](/oss/python/releases/langchain-v1)

{' ' * (self.list_depth - 1)}- [Migration guide](/oss/python/migrate/langchain-v1)

### [​](#langgraph)`langgraph`

{' ' * (self.list_depth - 1)}- [Release notes](/oss/python/releases/langgraph-v1)

{' ' * (self.list_depth - 1)}- [Migration guide](/oss/python/migrate/langgraph-v1)
If you encounter any issues or have feedback, please [open an issue](https://github.com/langchain-ai/docs/issues/new?template=01-langchain.yml) so we can improve. To view v0.x documentation, [go to the archived content](https://github.com/langchain-ai/langchain/tree/v0.3/docs/docs) and [API reference](https://reference.langchain.com/v0.3/python/).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/python/releases/changelog.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[VersioningPrevious](/oss/python/versioning)[What's new in LangChain v1Next](/oss/python/releases/langchain-v1)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
