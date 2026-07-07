# Unit testing

[Agent development](/oss/python/langchain/studio)[Test](/oss/python/langchain/test/index)

# Unit testing
Copy page

Test agent logic without API calls using fake chat models and in-memory persistence.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Unit tests exercise small, deterministic pieces of your agent in isolation. By replacing the real LLM with an in-memory fake (AKA fixture), you can script exact responses (text, tool calls, and errors) so tests are fast, free, and repeatable without API keys.

## [​](#mock-chat-model)Mock chat model

LangChain provides [`GenericFakeChatModel`](https://reference.langchain.com/python/langchain-core/language_models/fake_chat_models/GenericFakeChatModel) for mocking text responses. It accepts an iterator of responses ([`AIMessage`](https://reference.langchain.com/python/langchain-core/messages/ai/AIMessage) objects or strings) and returns one per invocation. It supports both regular and streaming usage.

```
from langchain_core.language_models.fake_chat_models import GenericFakeChatModel

model = GenericFakeChatModel(messages=iter([
 AIMessage(content="", tool_calls=[ToolCall(name="foo", args={"bar": "baz"}, id="call_1")]),
 "bar"
]))

model.invoke("hello")
# AIMessage(content='', ..., tool_calls=[{'name': 'foo', 'args': {'bar': 'baz'}, 'id': 'call_1', 'type': 'tool_call'}])

```

If we invoke the model again, it will return the next item in the iterator:

```
model.invoke("hello, again!")
# AIMessage(content='bar', ...)

```

## [​](#inmemorysaver-checkpointer)InMemorySaver checkpointer

To enable persistence during testing, you can use the [`InMemorySaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.memory.InMemorySaver) checkpointer. This allows you to simulate multiple turns to test state-dependent behavior:

```
from langgraph.checkpoint.memory import InMemorySaver

agent = create_agent(
 model,
 tools=[],
 checkpointer=InMemorySaver()
)

# First invocation
agent.invoke(
 {"messages": [HumanMessage(content="I live in Sydney, Australia")]},
 config={"configurable": {"thread_id": "session-1"}}
)

# Second invocation: the first message is persisted (Sydney location), so the model returns GMT+10 time
agent.invoke(
 {"messages": [HumanMessage(content="What's my local time?")]},
 config={"configurable": {"thread_id": "session-1"}}
)

```

## [​](#next-steps)Next steps

Learn how to test your agent with real model provider APIs in [Integration testing](/oss/python/langchain/test/integration-testing).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/test/unit-testing.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[TestPrevious](/oss/python/langchain/test)[Integration testingNext](/oss/python/langchain/test/integration-testing)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
