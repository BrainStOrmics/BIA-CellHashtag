# Quickstart

[Get started](/oss/python/langchain/install)

# Quickstart
Copy page

Build your first agent in minutesCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This quickstart shows you how to create a fully functional AI agent in just a few minutes.
**Using an AI coding assistant?**

{' ' * (self.list_depth - 1)}- Install the [LangChain Docs MCP server](/use-these-docs) to give your agent access to up-to-date LangChain documentation and examples.

{' ' * (self.list_depth - 1)}- Install [LangChain Skills](https://github.com/langchain-ai/langchain-skills) to improve your agent’s performance on LangChain ecosystem tasks.

## [​](#install-dependencies)Install dependencies

Install the following packages to follow along:
uvpipvenv
```
uv init
uv add langchain deepagents
uv sync

```

## [​](#set-up-api-keys)Set up API keys

Get an API key from [any supported model provider](/oss/python/integrations/providers/overview) (for example, Google Gemini or OpenAI).
Set the API keys, for example:

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Google Gemini
{' ' * (self.list_depth - 1)}- Claude (Anthropic)
{' ' * (self.list_depth - 1)}- OpenRouter
{' ' * (self.list_depth - 1)}- Fireworks
{' ' * (self.list_depth - 1)}- Baseten
{' ' * (self.list_depth - 1)}- Ollama
{' ' * (self.list_depth - 1)}- Azure
{' ' * (self.list_depth - 1)}- AWS Bedrock
{' ' * (self.list_depth - 1)}- HuggingFace
{' ' * (self.list_depth - 1)}- Other
```
export OPENAI_API_KEY="your-api-key"

```

```
export GOOGLE_API_KEY="your-api-key"

```

```
export ANTHROPIC_API_KEY="your-api-key"

```

```
export OPENROUTER_API_KEY="your-api-key"

```

```
export FIREWORKS_API_KEY="your-api-key"

```

```
export BASETEN_API_KEY="your-api-key"

```

```
# Local: Ollama must be running (https://ollama.com)
# Cloud: Set your Ollama API key for hosted inference
export OLLAMA_API_KEY="your-api-key"

```

```
export AZURE_OPENAI_API_KEY="your-api-key"
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com"
export AZURE_OPENAI_DEPLOYMENT_NAME="your-deployment"

```

```
export AWS_ACCESS_KEY_ID="your-access-key"
export AWS_SECRET_ACCESS_KEY="your-secret-key"
export AWS_REGION="us-east-1"

```

```
export HUGGINGFACEHUB_API_TOKEN="hf_..."

```
See the full list of supported [chat model integrations](/oss/python/integrations/chat).

## [​](#build-a-basic-agent)Build a basic agent

Start by creating a simple agent that can answer questions and call tools. The agent in this example uses the chosen language model, a basic weather function as a tool, and a simple prompt to guide its behavior:
OpenAIGoogle GeminiClaude (Anthropic)OpenRouterFireworksBasetenOllamaAzureAWS BedrockHuggingFace
```
from langchain.agents import create_agent

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

agent = create_agent(
 model="openai:gpt-5.4",
 tools=[get_weather],
 system_prompt="You are a helpful assistant",
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What's the weather in San Francisco?"}]}
)
print(result["messages"][-1].content_blocks)

```

When you run the code and prompt the agent to tell you about the weather in San Francisco, the agent uses that input and its available context.
The agent understands that you are asking about the weather for the city San Francisco and therefore calls the weather tool with the provided city name.
You can use [any supported model](/oss/python/integrations/providers/overview) by changing the model name and setting up the appropriate API key. Trace what’s happening inside your agent with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-quickstart)—follow the [tracing quickstart](/langsmith/trace-with-langchain) to get set up.

## [​](#build-a-real-world-agent)Build a real-world agent

In the following example you will build a research agent that can answer questions about text files.
Along the way you will explore the following concepts:

{' ' * (self.list_depth - 1)}- **Detailed system prompts** for better agent behavior

{' ' * (self.list_depth - 1)}- **Create tools** that integrate with external data

{' ' * (self.list_depth - 1)}- **Model configuration** for consistent responses

{' ' * (self.list_depth - 1)}- **Conversational memory** for chat-like interactions

{' ' * (self.list_depth - 1)}- **Deep Agents** for built-in features

{' ' * (self.list_depth - 1)}- **Testing** your agent

1[](#)

Define the system promptThe system prompt defines your agent’s role and behavior. Keep it specific and actionable:
```
SYSTEM_PROMPT = """You are a literary data assistant.

## Capabilities

- `fetch_text_from_url`: loads document text from a URL into the conversation.
Do not guess line counts or positions—ground them in tool results from the saved file."""

```
2[](#)

Create tools[Tools](/oss/python/langchain/tools) let a model interact with external systems by calling functions you define.
Tools can depend on [runtime context](/oss/python/langchain/runtime) and also interact with [agent memory](/oss/python/langchain/short-term-memory).This example uses a tool to load a document from a given URL:
```
import urllib.error
import urllib.request

from langchain.tools import tool

@tool
def fetch_text_from_url(url: str) -> str:
 """Fetch the document from a URL.
 """
 req = urllib.request.Request(
 url,
 headers={"User-Agent": "Mozilla/5.0 (compatible; quickstart-research/1.0)"},
 )
 try:
 with urllib.request.urlopen(req, timeout=120) as resp:
 raw = resp.read()
 except urllib.error.URLError as e:
 return f"Fetch failed: {e}"
 text = raw.decode("utf-8", errors="replace")
 return text

```
Tools should be well-documented: their name, description, and argument names become part of the model’s prompt.
LangChain’s [`@tool` decorator](https://reference.langchain.com/python/langchain-core/tools/convert/tool) adds metadata and enables runtime injection with the `ToolRuntime` parameter.
Learn more in the [tools guide](/oss/python/langchain/tools).3[](#)

Configure your modelSet up your [language model](/oss/python/langchain/models) with the right parameters for your use case. For example:OpenAIGoogle GeminiClaude (Anthropic)OpenRouterFireworksBasetenOllamaAzureAWS BedrockHuggingFace
```
from langchain.chat_models import init_chat_model

model = init_chat_model(
 "openai:gpt-5.4",
 temperature=0.5,
 timeout=300,
 max_tokens=25000,
)

```
Depending on the model and provider chosen, initialization parameters may vary; refer to their reference pages for details.4[](#)

Add memoryAdd [memory](/oss/python/langchain/short-term-memory) to your agent to maintain state across interactions. This allows
the agent to remember previous conversations and context.
```
from langgraph.checkpoint.memory import InMemorySaver

checkpointer = InMemorySaver()

```
In production, use a persistent checkpointer that saves message history to a database.
See [Add and manage memory](/oss/python/langgraph/add-memory#manage-short-term-memory) for more details.5[](#)

Create and run the agentNow assemble your agent with all the components and run it.There are two different frameworks for creating agents: LangChain agents and deep agents.
Both LangChain and deep agents provide you with fine-grained control over tools, memory, and more.
The main difference between both is that deep agents come with a range of commonly useful capabilities already built in, such as planning, file system tools, and subagents.Use deep agents when you want maximum capability with minimal setup; choose LangChain agents when you need fine-grained control.Since the code invokes the model with the entire text from The Great Gatsby, it uses a large amount of tokens.You can view example output in the next step.Let’s try both:
```
from langchain.agents import create_agent
from deepagents import create_deep_agent

agent = create_agent(
 model=model,
 tools=[fetch_text_from_url],
 system_prompt=SYSTEM_PROMPT,
 checkpointer=checkpointer,
)

deep_agent = create_deep_agent(
 model=model,
 tools=[fetch_text_from_url],
 system_prompt=SYSTEM_PROMPT,
 checkpointer=checkpointer,
)

content = f"""Project Gutenberg hosts a full plain-text copy of F. Scott Fitzgerald's The Great Gatsby.
URL: https://www.gutenberg.org/files/64317/64317-0.txt

Answer as much as you can:

1) How many lines in the complete Gutenberg file contain the substring `Gatsby` (count lines, not occurrences within a line, each line ends with a line break).
2) The 1-based line number of the first line in the file that contains `Daisy`.
3) A two-sentence neutral synopsis.

Do your best on (1) and (2). If at any point you realize you cannot **verify** an exact answer with
your available tools and reasoning, do not fabricate numbers: use `null` for that field and spell out
the limitation in `how_you_computed_counts`. If you encounter any errors please report what the error was and what the error message was."""

agent_result = agent.invoke(
 {"messages": [{"role": "user", "content": content}]},
 config={"configurable": {"thread_id": "great-gatsby-lc"}},
)
deep_agent_result = deep_agent.invoke(
 {"messages": [{"role": "user", "content": content}]},
 config={"configurable": {"thread_id": "great-gatsby-da"}},
)
print(agent_result["messages"][-1].content_blocks)
print("\n")
print(deep_agent_result["messages"][-1].content_blocks)

```

Show Full example code
```
import urllib.error
import urllib.request

from langchain.agents import create_agent
from deepagents import create_deep_agent
from langchain.chat_models import init_chat_model
from langchain.tools import tool
from langgraph.checkpoint.memory import InMemorySaver

SYSTEM_PROMPT = """You are a literary data assistant.

## Capabilities

- `fetch_text_from_url`: loads document text from a URL into the conversation.
Do not guess line counts or positions—ground them in tool results from the saved file."""

@tool
def fetch_text_from_url(url: str) -> str:
 """Fetch the document from a URL.
 """
 req = urllib.request.Request(
 url,
 headers={"User-Agent": "Mozilla/5.0 (compatible; quickstart-research/1.0)"},
 )
 try:
 with urllib.request.urlopen(req, timeout=120) as resp:
 raw = resp.read()
 except urllib.error.URLError as e:
 return f"Fetch failed: {e}"
 text = raw.decode("utf-8", errors="replace")
 return text

model = init_chat_model(
 "gemini-3.1-pro-preview",
 model_provider="google-genai",
 temperature=0.5,
 timeout=600,
 max_tokens=25000,
 streaming=True,
)

checkpointer = InMemorySaver()

agent = create_agent(
 model=model,
 tools=[fetch_text_from_url],
 system_prompt=SYSTEM_PROMPT,
 checkpointer=checkpointer,
)

deep_agent = create_deep_agent(
 model=model,
 tools=[fetch_text_from_url],
 system_prompt=SYSTEM_PROMPT,
 checkpointer=checkpointer,
)

content = f"""Project Gutenberg hosts a full plain-text copy of F. Scott Fitzgerald's The Great Gatsby.
URL: https://www.gutenberg.org/files/64317/64317-0.txt

Answer as much as you can:

1) How many lines in the complete Gutenberg file contain the substring `Gatsby` (count lines, not occurrences within a line, each line ends with a line break).
2) The 1-based line number of the first line in the file that contains `Daisy`.
3) A two-sentence neutral synopsis.

Do your best on (1) and (2). If at any point you realize you cannot **verify** an exact answer with
your available tools and reasoning, do not fabricate numbers: use `null` for that field and spell out
the limitation in `how_you_computed_counts`. If you encounter any errors please report what the error was and what the error message was."""

agent_result = agent.invoke(
 {"messages": [{"role": "user", "content": content}]},
 config={"configurable": {"thread_id": "great-gatsby-lc"}},
)
deep_agent_result = deep_agent.invoke(
 {"messages": [{"role": "user", "content": content}]},
 config={"configurable": {"thread_id": "great-gatsby-da"}},
)
print(agent_result["messages"][-1].content_blocks)
print("\n")
print(deep_agent_result["messages"][-1].content_blocks)

```
6[](#)

Review the resultsThe results will differ based on the model and the execution.
{' ' * (self.list_depth - 1)}- LangChain agents
{' ' * (self.list_depth - 1)}- Deep agents
```
**1) Number of lines containing `Gatsby`:** `null`

**2) First line containing `Daisy`:** `null`

**3) Synopsis:**
The Great Gatsby follows the mysterious millionaire Jay Gatsby and his obsession with reuniting with his former lover, Daisy Buchanan, as narrated by his neighbor Nick Carraway. Set against the backdrop of the Roaring Twenties on Long Island, the novel explores themes of wealth, class, and the elusive nature of the American Dream.

**how_you_computed_counts:**
I successfully fetched the full text of the eBook using the `fetch_text_from_url` tool. However, because I do not have access to a code execution environment (like Python) or text-processing tools (like `grep`), I cannot deterministically split the text by line breaks, iterate through the thousands of lines, and verify the exact line numbers or match counts. LLMs cannot reliably perform exact line-counting or indexing over massive texts within their context window without external computational tools. As instructed, rather than fabricating or guessing a number, I have output `null` for the exact counts and positions.

```
See all 10 lines
```
Based on the text fetched directly from the Gutenberg URL and analyzed using filesystem search tools, here are the answers to your questions:

**1) Lines containing the substring `Gatsby`**
**258** lines contain the exact substring `Gatsby`.

**2) First line containing `Daisy`**
Line **181** is the first line in the file that contains the exact substring `Daisy`.
*(For context, the line reads: "Buchanans. Daisy was my second cousin once removed, and I’d known Tom")*

**3) Two-sentence neutral synopsis**
*The Great Gatsby* follows the mysterious millionaire Jay Gatsby and his obsessive pursuit to reunite with his former lover, Daisy Buchanan, in 1920s Long Island. The story is narrated by Nick Carraway, who observes the tragic consequences of Gatsby's relentless ambition and the shallow materialism of the era's wealthy elite.

***

**How counts were computed:**
When fetching the document from the URL, the file was too large for the standard output and was automatically saved to the local filesystem by the system (`/large_tool_results/x246ax2x`). I then used the `grep` tool to search the saved file for the exact literal substrings `Gatsby` and `Daisy`. The `grep` tool returned every matching line along with its 1-based line number. I manually counted the exact number of lines returned for `Gatsby` (which totaled 258) and identified the first line number returned for `Daisy` (which was 181). I also verified there were no uppercase variations (`GATSBY` or `DAISY`) that would have been missed. No errors were encountered during this process.

```
See all 17 linesIf you look at the output on both tabs, you notice that the LangChain agent provided answers but they are estimates. The agent lacks the tools to answer this question. You may also get errors that the prompt is too long.The deep agent, on the other hand can:

{' ' * (self.list_depth - 1)}- **Plans its approach** using the built-in [`write_todos`](/oss/python/deepagents/harness#task-planning) tool to break down the research task.

{' ' * (self.list_depth - 1)}- **Loads the file** by calling the `fetch_text_from_url` tool to gather information.

{' ' * (self.list_depth - 1)}- **Manages context** by using file system tools ([`grep`](/oss/python/deepagents/harness#virtual-filesystem-access) and [`read_file`](/oss/python/deepagents/harness#virtual-filesystem-access)).

{' ' * (self.list_depth - 1)}- **Spawns subagents** as needed to delegate complex subtasks to specialized subagents.
For LangChain agents, you must implement more capabilities to get a similar level of service and can customize them along the way as needed.

## [​](#trace-agent-calls)Trace agent calls

Most interesting applications you build with LangChain make many calls to LLMs. As these applications get more complex, it becomes important to be able to inspect what exactly is going on inside your agent. The best way to do this is with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-quickstart).
Sign up for a [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-quickstart) account and set these to start logging traces:

```
export LANGSMITH_TRACING="true"
export LANGSMITH_API_KEY="..."

```

Once set, run your script again and then inspect what happened during your agent calls on [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-quickstart) .
To learn more about tracing your agent with LangSmith, see the [LangSmith documentation](/langsmith/trace-with-langchain).

## [​](#next-steps)Next steps

You now have agents that can:

{' ' * (self.list_depth - 1)}- **Understand context** and remember conversations

{' ' * (self.list_depth - 1)}- **Use tools** intelligently

{' ' * (self.list_depth - 1)}- **Provide structured responses** in a consistent format

{' ' * (self.list_depth - 1)}- **Handle user-specific information** through context

{' ' * (self.list_depth - 1)}- **Maintain conversation state** across interactions

{' ' * (self.list_depth - 1)}- **Plan, research, and synthesize** (deep agents only)

Continue with:

{' ' * (self.list_depth - 1)}- **LangChain agents**: [Add and manage memory](/oss/python/langgraph/add-memory#manage-short-term-memory), [deploy to production](/oss/python/langgraph/deploy)

{' ' * (self.list_depth - 1)}- **Deep Agents**: [Customization options](/oss/python/deepagents/customization), [persistent memory](/oss/python/deepagents/long-term-memory), [deploy to production](/oss/python/langgraph/deploy)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/quickstart.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Install LangChainPrevious](/oss/python/langchain/install)[ChangelogNext](/oss/python/langchain/changelog-py)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
