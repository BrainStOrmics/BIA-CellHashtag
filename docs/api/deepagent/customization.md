# Customize Deep Agents

[Get started](/oss/python/deepagents/quickstart)

# Customize Deep Agents
Copy page

Learn how to customize Deep Agents with system prompts, tools, subagents, and moreCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Build the harness around your goal. `create_deep_agent` gives you a production-ready foundation: connect it to your data, shape its behavior, and add the capabilities your use case needs.

```
from deepagents import create_deep_agent

agent = create_deep_agent(
 model="anthropic:claude-sonnet-4-6",
 system_prompt="You are a helpful assistant.",
 tools=[search, fetch_url],
 memory=["./AGENTS.md"],
 skills=["./skills/"],
)

```

| 
 | Parameter | What it does
 | [`model=`](#model) | Which model to use
 | [`system_prompt=`](#system-prompt) | Custom instructions for the agent
 | [`tools=`](#tools) | Domain tools the agent can call
 | [`memory=`](#memory) | AGENTS.md files loaded at startup
 | [`skills=`](#skills) | Skills directory for on-demand knowledge
 | [`backend=`](#backends) | Filesystem backend (StateBackend by default)
 | [`permissions=`](/oss/python/deepagents/permissions) | Path-level access control for the filesystem
 | [`subagents=`](#subagents) | Custom subagents for delegated tasks
 | [`middleware=`](#middleware) | Extra middleware to add to the stack
 | [`interrupt_on=`](#human-in-the-loop) | Pause before tool calls for human approval
 | [`response_format=`](#structured-output) | Structured output schema
 | [profiles](#profiles) | Per-model defaults as a reusable bundle

Full function signature
```
create_deep_agent(
 model: str | BaseChatModel | None = None,
 tools: Sequence[BaseTool | Callable | dict[str, Any]] | None = None,
 *,
 system_prompt: str | SystemMessage | None = None,
 middleware: Sequence[AgentMiddleware] = (),
 subagents: Sequence[SubAgent | CompiledSubAgent | AsyncSubAgent] | None = None,
 skills: list[str] | None = None,
 memory: list[str] | None = None,
 permissions: list[FilesystemPermission] | None = None,
 backend: BackendProtocol | BackendFactory | None = None,
 interrupt_on: dict[str, bool | InterruptOnConfig] | None = None,
 response_format: ResponseFormat[ResponseT] | type[ResponseT] | dict[str, Any] | None = None,
 context_schema: type[ContextT] | None = None,
 checkpointer: Checkpointer | None = None,
 store: BaseStore | None = None,
 debug: bool = False,
 name: str | None = None,
 cache: BaseCache | None = None
) -> CompiledStateGraph[AgentState[ResponseT], ContextT, _InputAgentState, _OutputAgentState[ResponseT]]

```

For the full parameter list, see the [`create_deep_agent`](https://reference.langchain.com/python/deepagents/graph/create_deep_agent) API reference. To compose a fully custom harness from scratch, see [Configure the harness](/oss/python/langchain/agents#configure-the-harness) or follow the step-by-step [Build a deep agent from scratch](/oss/python/langchain/deep-agent-from-scratch) guide.

## [​](#model)Model

Pass a `model` string in `provider:model` format, or an initialized model instance. See [supported models](/oss/python/deepagents/models#supported-models) for all providers and [suggested models](/oss/python/deepagents/models#suggested-models) for tested recommendations.
Use the `provider:model` format (for example `openai:gpt-5.4`) to quickly switch between models.

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- Azure
{' ' * (self.list_depth - 1)}- Google Gemini
{' ' * (self.list_depth - 1)}- AWS Bedrock
{' ' * (self.list_depth - 1)}- HuggingFace
{' ' * (self.list_depth - 1)}- Other👉 Read the [OpenAI chat model integration docs](/oss/python/integrations/chat/openai)
```
pip install -U "langchain[openai]"

```
default parametersinit_chat_modelModel Class
```
import os
from deepagents import create_deep_agent

os.environ["OPENAI_API_KEY"] = "sk-..."

agent = create_deep_agent(model="openai:gpt-5.4")
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
👉 Read the [Anthropic chat model integration docs](/oss/python/integrations/chat/anthropic)
```
pip install -U "langchain[anthropic]"

```
default parametersinit_chat_modelModel Class
```
import os
from deepagents import create_deep_agent

os.environ["ANTHROPIC_API_KEY"] = "sk-..."

agent = create_deep_agent(model="anthropic:claude-sonnet-4-6")
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
👉 Read the [Azure chat model integration docs](/oss/python/integrations/chat/azure_chat_openai)
```
pip install -U "langchain[openai]"

```
default parametersinit_chat_modelModel Class
```
import os
from deepagents import create_deep_agent

os.environ["AZURE_OPENAI_API_KEY"] = "..."
os.environ["AZURE_OPENAI_ENDPOINT"] = "..."
os.environ["OPENAI_API_VERSION"] = "2025-03-01-preview"

agent = create_deep_agent(model="azure_openai:gpt-5.4")
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
👉 Read the [Google GenAI chat model integration docs](/oss/python/integrations/chat/google_generative_ai)
```
pip install -U "langchain[google-genai]"

```
default parametersinit_chat_modelModel Class
```
import os
from deepagents import create_deep_agent

os.environ["GOOGLE_API_KEY"] = "..."

agent = create_deep_agent(model="google_genai:gemini-3.5-flash")
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
👉 Read the [AWS Bedrock chat model integration docs](/oss/python/integrations/chat/bedrock)
```
pip install -U "langchain[aws]"

```
default parametersinit_chat_modelModel Class
```
from deepagents import create_deep_agent

# Follow the steps here to configure your credentials:
# https://docs.aws.amazon.com/bedrock/latest/userguide/getting-started.html

agent = create_deep_agent(
 model="anthropic.claude-sonnet-4-6",
 model_provider="bedrock_converse",
)
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
👉 Read the [HuggingFace chat model integration docs](/oss/python/integrations/chat/huggingface)
```
pip install -U "langchain[huggingface]"

```
default parametersinit_chat_modelModel Class
```
import os
from deepagents import create_deep_agent

os.environ["HUGGINGFACEHUB_API_TOKEN"] = "hf_..."

agent = create_deep_agent(
 model="microsoft/Phi-3-mini-4k-instruct",
 model_provider="huggingface",
 temperature=0.7,
 max_tokens=1024,
)
# this calls init_chat_model for the specified model with default parameters
# to use specific model parameters, use init_chat_model directly

```
Pass any [supported model string](/oss/python/deepagents/models#supported-models), or an initialized model instance:model stringinit_chat_modelmodel class
```
from deepagents import create_deep_agent

agent = create_deep_agent(model="provider:model-name")

```

Chat models automatically retry transient API failures (with exponential backoff). For defaults, limits, and code samples for tuning `max_retries` / `timeout` live on the LangChain [Models](/oss/python/langchain/models#connection-resilience) page.

## [​](#tools)Tools

In addition to [built-in tools](/oss/python/deepagents/overview#core-capabilities) for planning, file management, and subagent spawning, you can provide custom tools:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
import os
from typing import Literal
from tavily import TavilyClient
from deepagents import create_deep_agent

tavily_client = TavilyClient(api_key=os.environ["TAVILY_API_KEY"])

def internet_search(
 query: str,
 max_results: int = 5,
 topic: Literal["general", "news", "finance"] = "general",
 include_raw_content: bool = False,
):
 """Run a web search"""
 return tavily_client.search(
 query,
 max_results=max_results,
 include_raw_content=include_raw_content,
 topic=topic,
 )

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[internet_search],
)

```

## [​](#system-prompt)System prompt

Deep Agents come with a built-in system prompt. A deep agent’s value comes from the orchestration layer the SDK provides on top of the model—planning, virtual-filesystem tools, and subagents—and the model needs to know those exist and when to reach for them. The built-in prompt teaches the agent how to use that scaffolding so you don’t have to re-derive it for every project; tweak it through a [profile](/oss/python/deepagents/profiles#harness-profiles) or your own `system_prompt=` rather than copying it verbatim.
When middleware add special tools, like the filesystem tools, it appends them to the system prompt.
Each deep agent should also include a custom system prompt specific to its specific use case:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from deepagents import create_deep_agent

research_instructions = """\
You are an expert researcher. Your job is to conduct \
thorough research, and then write a polished report. \
"""

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 system_prompt=research_instructions,
)

```

### [​](#prompt-assembly)Prompt assembly

Deep Agents builds the system prompt from up to four named parts so that caller-supplied instructions, the SDK’s built-in agent guidance, and any model-specific [profile](/oss/python/deepagents/profiles) overrides can coexist with predictable precedence. Without this layering, a profile suffix tuned for Claude (for example) could overwrite or be overwritten by your `system_prompt=` argument depending on call order; the named slots make the ordering explicit and stable.
In practice, most callers only encounter two slots: `USER` (your `system_prompt=`) and `BASE` (the SDK default). Selecting a model with a built-in profile—Anthropic or OpenAI today—adds a `SUFFIX`. The full four-part assembly is mainly relevant when you author a custom `HarnessProfile` or debug why a profile’s text appears where it does.
The four named parts (each may be absent):

| 
 | Name | Source | Notes
 | `USER` | `system_prompt=` argument to `create_deep_agent` | `str` or `SystemMessage`; omitted when unset.
 | `BASE` | The SDK default (`BASE_AGENT_PROMPT`) | Always present unless replaced by a profile’s `CUSTOM`.
 | `CUSTOM` | [`HarnessProfile.base_system_prompt`](/oss/python/deepagents/profiles#harness-profiles) | Replaces `BASE` outright when a matching profile sets it.
 | `SUFFIX` | [`HarnessProfile.system_prompt_suffix`](/oss/python/deepagents/profiles#harness-profiles) | Appended last when a matching profile sets it.
The order is always **`USER` -> (`BASE` or `CUSTOM`) -> `SUFFIX`**, joined by blank lines (`\n\n`). Two invariants follow:

{' ' * (self.list_depth - 1)}- **`USER` is always at the front.** The caller’s text precedes any SDK or profile content, so persona/instructions take precedence regardless of which model is selected.

{' ' * (self.list_depth - 1)}- **`SUFFIX` is always at the end.** Profile suffixes sit closest to the conversation history, where model-tuning guidance lands most reliably.

Assembled shapes (✓ = field is set, - = field is unset):

| 
 | `system_prompt=` | profile `base_system_prompt` (`CUSTOM`) | profile `system_prompt_suffix` (`SUFFIX`) | Final assembled system prompt
 | `None` | - | - | `BASE`
 | `None` | - | ✓ | `BASE` + `SUFFIX`
 | `None` | ✓ | - | `CUSTOM`
 | `None` | ✓ | ✓ | `CUSTOM` + `SUFFIX`
 | `str` | - | - | `USER` + `BASE`
 | `str` | - | ✓ | `USER` + `BASE` + `SUFFIX`
 | `str` | ✓ | - | `USER` + `CUSTOM`
 | `str` | ✓ | ✓ | `USER` + `CUSTOM` + `SUFFIX`
Worked example—built-in profiles (Anthropic, OpenAI) ship only a `system_prompt_suffix`, so a typical call lands in the `str` + `-` + `✓` row:

```
agent = create_deep_agent(
 model="anthropic:claude-sonnet-4-6",
 system_prompt="You are a customer-support agent for ACME Corp.",
)
# Final = USER + BASE + SUFFIX
# = "You are a customer-support agent for ACME Corp."
# + "\n\n"
# + BASE_AGENT_PROMPT
# + "\n\n"
# + <Claude-specific guidance>

```

Passing a `SystemMessage` (rather than a string) triggers a different concatenation path: the right-hand assembly (`BASE`-or-`CUSTOM` plus any `SUFFIX`) is appended as an additional text content block onto the message’s existing `content_blocks`. The same logical ordering applies (caller blocks first), and any `cache_control` markers on the caller’s blocks are preserved—useful for placing explicit Anthropic prompt-cache breakpoints.

Subagent promptsThe same overlay rules apply to declarative [subagents](/oss/python/deepagents/subagents)—each subagent re-runs profile resolution against **its own model**, then applies the resolved profile’s `base_system_prompt` / `system_prompt_suffix` to its authored `system_prompt`. The subagent’s `system_prompt` plays the `BASE` role; `CUSTOM` and `SUFFIX` come from the profile that matches the subagent’s model (which may differ from the main agent’s profile).

| 
 | `spec["system_prompt"]` | profile `base_system_prompt` (`CUSTOM`) | profile `system_prompt_suffix` (`SUFFIX`) | Final subagent system prompt
 | authored | - | - | authored
 | authored | - | ✓ | authored + `SUFFIX`
 | authored | ✓ | - | `CUSTOM`
 | authored | ✓ | ✓ | `CUSTOM` + `SUFFIX`There is no `USER` segment for subagents—the spec’s authored `system_prompt` is the closest analog and stays in the `BASE` slot. A profile that ships only a `system_prompt_suffix` (the common case for built-in Anthropic / OpenAI profiles) just appends to whatever the subagent author wrote; a profile that sets `base_system_prompt` will *replace* the authored prompt outright, so reach for that field deliberately.

General-purpose subagent promptThe auto-added [general-purpose subagent](/oss/python/deepagents/subagents#the-general-purpose-subagent) follows the same overlay rules with one extra layer: the GP base prompt is resolved as **`general_purpose_subagent.system_prompt` (if set) -> `HarnessProfile.base_system_prompt` (if set) -> SDK GP default**. The profile suffix layers on top either way.The two override fields can both carry a base-prompt replacement, but they are not interchangeable. `general_purpose_subagent.system_prompt` is GP-specific configuration; `base_system_prompt` is a global override that primarily targets the main agent. When both are set, the **GP-specific intent wins for the GP subagent** so a user tuning both fields never sees their GP override silently dropped:
```
register_harness_profile(
 "anthropic",
 HarnessProfile(
 base_system_prompt="You are ACME's support orchestrator.", # main agent
 general_purpose_subagent=GeneralPurposeSubagentProfile(
 system_prompt="You are a research subagent. Cite sources.", # GP subagent
 ),
 system_prompt_suffix="Always think step by step.",
 ),
)

```

| 
 | Stack | Final system prompt
 | Main agent | `"You are ACME's support orchestrator." + SUFFIX`
 | GP subagent | `"You are a research subagent. Cite sources." + SUFFIX`If `general_purpose_subagent.system_prompt` is unset, the GP subagent falls back to `base_system_prompt` (when set) and finally to the SDK GP default.

## [​](#middleware)Middleware

Deep Agents support any [middleware](/oss/python/langchain/middleware/overview), including the built-in middleware listed below, prebuilt middleware from LangChain, provider-specific middleware, and custom middleware you write yourself. Pass middleware to the `middleware` argument of `create_deep_agent`.
By default, Deep Agents have access to the following middleware:

{' ' * (self.list_depth - 1)}- [`TodoListMiddleware`](https://reference.langchain.com/python/langchain/agents/middleware/todo/TodoListMiddleware): Tracks and manages todo lists for organizing agent tasks and work

{' ' * (self.list_depth - 1)}- [`FilesystemMiddleware`](https://reference.langchain.com/python/deepagents/middleware/filesystem/FilesystemMiddleware): Handles file system operations such as reading, writing, and navigating directories

{' ' * (self.list_depth - 1)}- [`SubAgentMiddleware`](https://reference.langchain.com/python/deepagents/middleware/subagents/SubAgentMiddleware): Spawns and coordinates subagents for delegating tasks to specialized agents

{' ' * (self.list_depth - 1)}- [`SummarizationMiddleware`](https://reference.langchain.com/python/langchain/agents/middleware/summarization/SummarizationMiddleware): Condenses message history to stay within context limits when conversations grow long

{' ' * (self.list_depth - 1)}- [`AnthropicPromptCachingMiddleware`](https://reference.langchain.com/python/langchain-anthropic/middleware/prompt_caching/AnthropicPromptCachingMiddleware): Automatic reduction of redundant token processing when using Anthropic models

{' ' * (self.list_depth - 1)}- [`PatchToolCallsMiddleware`](https://reference.langchain.com/python/deepagents/middleware/patch_tool_calls/PatchToolCallsMiddleware): Automatic message history fixes when tool calls are interrupted or cancelled before receiving results

If you are using memory, skills, or human-in-the-loop, the following middleware is also included:

{' ' * (self.list_depth - 1)}- [`MemoryMiddleware`](https://reference.langchain.com/python/deepagents/middleware/memory/MemoryMiddleware): Persists and retrieves conversation context across sessions when the `memory` argument is provided

{' ' * (self.list_depth - 1)}- [`SkillsMiddleware`](https://reference.langchain.com/python/deepagents/middleware/skills/SkillsMiddleware): Enables custom skills when the `skills` argument is provided

{' ' * (self.list_depth - 1)}- `HumanInTheLoopMiddleware`: Pauses for human approval or input at specified points when the `interruptOn` argument is provided

### [​](#prebuilt-middleware)Prebuilt middleware

LangChain exposes additional prebuilt middleware that let you add-on various features, such as retries, fallbacks, or PII detection. See [Prebuilt middleware](/oss/python/langchain/middleware/built-in) for more.
The `deepagents` library also exposes [`create_summarization_tool_middleware`](https://reference.langchain.com/python/deepagents/middleware/summarization/create_summarization_tool_middleware), enabling agents to trigger summarization at opportune times—such as between tasks—instead of at fixed token intervals. For more detail, see [Summarization](/oss/python/deepagents/context-engineering#summarization).

### [​](#provider-specific-middleware)Provider-specific middleware

For provider-specific middleware that is optimized for specific LLM providers, see [Official integrations](/oss/python/integrations/middleware#official-integrations) and [Community integrations](/oss/python/integrations/middleware#community-integrations).

### [​](#custom-middleware)Custom middleware

You can provide additional middleware to extend functionality, add tools, or implement custom hooks:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from langchain.agents.middleware import wrap_tool_call
from langchain.tools import tool
from deepagents import create_deep_agent

@tool
def get_weather(city: str) -> str:
 """Get the weather in a city."""
 return f"The weather in {city} is sunny."

call_count = [0] # Use list to allow modification in nested function

@wrap_tool_call
def log_tool_calls(request, handler):
 """Intercept and log every tool call - demonstrates cross-cutting concern."""
 call_count[0] += 1
 tool_name = request.name if hasattr(request, "name") else str(request)

 print(f"[Middleware] Tool call #{call_count[0]}: {tool_name}")
 print(f"[Middleware] Arguments: {request.args if hasattr(request, 'args') else 'N/A'}")

 # Execute the tool call
 result = handler(request)

 # Log the result
 print(f"[Middleware] Tool call #{call_count[0]} completed")

 return result

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[get_weather],
 middleware=[log_tool_calls],
)

```

**Do not mutate attributes after initialization**If you need to track values across hook invocations (for example, counters or accumulated data), use graph state.
Graph state is scoped to a thread by design, so updates are safe under concurrency.**Do this:**
```
from langchain.agents.middleware import AgentMiddleware

class CustomMiddleware(AgentMiddleware):
 def __init__(self):
 pass

 def before_agent(self, state, runtime):
 return {"x": state.get("x", 0) + 1} # Update graph state instead

```
Do **not** do this:
```
class CustomMiddlewareBad(AgentMiddleware):
 def __init__(self):
 self.x = 1

 def before_agent(self, state, runtime):
 self.x += 1 # Mutation causes race conditions

```
Mutation in place, such as modifying `self.x` in `before_agent` or changing other shared values in hooks, can lead to subtle bugs and race conditions because many operations run concurrently (subagents, parallel tools, and parallel invocations on different threads).For full details on extending state with custom properties, see [Custom middleware - Custom state schema](/oss/python/langchain/middleware/custom#custom-state-schema).
If you must use mutation in custom middleware, consider what happens when subagents, parallel tools, or concurrent agent invocations run at the same time.

### [​](#interpreters)Interpreters

Use [interpreters](/oss/python/deepagents/interpreters) to add an `eval` tool that runs JavaScript in a scoped QuickJS runtime. Interpreters are useful when the agent needs to compose tools programmatically, batch work, handle errors in code, or transform structured data without a full shell environment.
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from deepagents import create_deep_agent
from langchain_quickjs import CodeInterpreterMiddleware

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 middleware=[CodeInterpreterMiddleware()],
)

```

For setup, programmatic tool calling, interpreter skills, and limits, see [Interpreters](/oss/python/deepagents/interpreters).

## [​](#subagents)Subagents

To isolate detailed work and avoid context bloat, use subagents:

```
import os
from typing import Literal

from deepagents import create_deep_agent
from tavily import TavilyClient

tavily_client = TavilyClient(api_key=os.environ["TAVILY_API_KEY"])

def internet_search(
 query: str,
 max_results: int = 5,
 topic: Literal["general", "news", "finance"] = "general",
 include_raw_content: bool = False,
):
 """Run a web search"""
 return tavily_client.search(
 query,
 max_results=max_results,
 include_raw_content=include_raw_content,
 topic=topic,
 )

research_subagent = {
 "name": "research-agent",
 "description": "Used to research more in depth questions",
 "system_prompt": "You are a great researcher",
 "tools": [internet_search],
 "model": "openai:gpt-5.4", # Optional override, defaults to main agent model
}
subagents = [research_subagent]

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 subagents=subagents,
)

```

For more information, see [Subagents](/oss/python/deepagents/subagents).

## [​](#backends)Backends

Tools for a deep agent can make use of virtual file systems to store, access, and edit files. By default, deep agents use a [`StateBackend`](https://reference.langchain.com/python/deepagents/backends/state/StateBackend).
If you are using [skills](#skills) or [memory](#memory), you must add the expected skill or memory files to the backend before creating the agent.

{' ' * (self.list_depth - 1)}- StateBackend
{' ' * (self.list_depth - 1)}- FilesystemBackend
{' ' * (self.list_depth - 1)}- LocalShellBackend
{' ' * (self.list_depth - 1)}- StoreBackend
{' ' * (self.list_depth - 1)}- ContextHubBackend
{' ' * (self.list_depth - 1)}- CompositeBackendA thread-scoped filesystem backend stored in `langgraph` state.Files persist across turns within a thread (via your checkpointer) and are not shared across threads.
```
from deepagents import create_deep_agent
from deepagents.backends import StateBackend

# By default we provide a StateBackend
agent = create_deep_agent(model="google_genai:gemini-3.5-flash")

# Under the hood, it looks like
agent2 = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=StateBackend(),
)

```
The local machine’s filesystem.This backend grants agents direct filesystem read/write access.
Use with caution and only in appropriate environments.
For more information, see [`FilesystemBackend`](/oss/python/deepagents/backends#filesystembackend-local-disk).
```
from deepagents import create_deep_agent
from deepagents.backends import FilesystemBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=FilesystemBackend(root_dir=".", virtual_mode=True),
)

```
Wrap `FilesystemBackend` in a `CompositeBackend` to prevent internal agent data (offloaded tool results, conversation history) from being written to disk alongside your project files. See the [recommended pattern](/oss/python/deepagents/backends#filesystembackend-local-disk).A filesystem with shell execution directly on the host. Provides filesystem tools plus the `execute` tool for running commands.This backend grants agents direct filesystem read/write access **and** unrestricted shell execution on your host.
Use with extreme caution and only in appropriate environments.
For more information, see [`LocalShellBackend`](/oss/python/deepagents/backends#localshellbackend-local-shell).
```
from deepagents import create_deep_agent
from deepagents.backends import LocalShellBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=LocalShellBackend(root_dir=".", virtual_mode=True, env={"PATH": "/usr/bin:/bin"}),
)

```
A filesystem that provides long-term storage that is *persisted across threads*.
```
from deepagents import create_deep_agent
from deepagents.backends import StoreBackend
from langgraph.store.memory import InMemoryStore

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
 ),
 store=InMemoryStore(), # Good for local dev; omit for LangSmith Deployment
)

```
When deploying to [LangSmith Deployment](/langsmith/deployment), omit the `store` parameter. The platform automatically provisions a store for your agent.The `namespace` parameter controls data isolation. For multi-user deployments, always set a [namespace factory](/oss/python/deepagents/backends#namespace-factories) to isolate data per user or tenant.Durable filesystem storage in a LangSmith Hub repo.
```
from deepagents import create_deep_agent
from deepagents.backends import ContextHubBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=ContextHubBackend("my-agent"),
)

```
For more details, see [`ContextHubBackend`](/oss/python/deepagents/backends#contexthubbackend).A flexible backend where you can specify different routes in the filesystem to point towards different backends.
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend
from langgraph.store.memory import InMemoryStore

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(namespace=lambda _rt: ("memories",)),
 },
 ),
 store=InMemoryStore(), # Store passed to create_deep_agent, not backend
)

```

For more information, see [Backends](/oss/python/deepagents/backends).

### [​](#sandboxes)Sandboxes

Sandboxes are specialized [backends](/oss/python/deepagents/backends) that run agent code in an isolated environment with their own filesystem and an `execute` tool for shell commands.
Use a sandbox backend when you want your deep agent to write files, install dependencies, and run commands without changing anything on your local machine.
You configure sandboxes by passing a sandbox backend to `backend` when creating your deep agent:

{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- LangSmithpipuv
```
pip install langchain-modal

```

```
import modal
from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_modal import ModalSandbox

app = modal.App.lookup("your-app")
modal_sandbox = modal.Sandbox.create(app=app)
backend = ModalSandbox(sandbox=modal_sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)
try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 modal_sandbox.terminate()

```
pipuv
```
pip install langchain-runloop

```

```
import os

from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_runloop import RunloopSandbox
from runloop_api_client import RunloopSDK

client = RunloopSDK(bearer_token=os.environ["RUNLOOP_API_KEY"])

devbox = client.devbox.create()
backend = RunloopSandbox(devbox=devbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)

try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 devbox.shutdown()

```
pipuv
```
pip install langchain-daytona

```

```
from daytona import Daytona
from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_daytona import DaytonaSandbox

sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)

try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 sandbox.stop()

```
LangSmith sandboxes are currently in private beta.pipuv
```
pip install "langsmith[sandbox]"

```

```
from deepagents import create_deep_agent
from deepagents.backends import LangSmithSandbox
from langchain_anthropic import ChatAnthropic
from langsmith.sandbox import SandboxClient

client = SandboxClient()
ls_sandbox = client.create_sandbox()
backend = LangSmithSandbox(sandbox=ls_sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)
try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 client.delete_sandbox(ls_sandbox.name)

```

For more information, see [Sandboxes](/oss/python/deepagents/sandboxes).

## [​](#human-in-the-loop)Human-in-the-loop

Some tool operations may be sensitive and require human approval before execution.
You can configure the approval for each tool:

```
from langchain.tools import tool
from deepagents import create_deep_agent
from langgraph.checkpoint.memory import MemorySaver

@tool
def remove_file(path: str) -> str:
 """Delete a file from the filesystem."""
 return f"Deleted {path}"

@tool
def fetch_file(path: str) -> str:
 """Read a file from the filesystem."""
 return f"Contents of {path}"

@tool
def notify_email(to: str, subject: str, body: str) -> str:
 """Send an email."""
 return f"Sent email to {to}"

# Checkpointer is REQUIRED for human-in-the-loop
checkpointer = MemorySaver()

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[remove_file, fetch_file, notify_email],
 interrupt_on={
 "remove_file": True, # Default: approve, edit, reject, respond
 "fetch_file": False, # No interrupts needed
 "notify_email": {"allowed_decisions": ["approve", "reject"]}, # No editing
 },
 checkpointer=checkpointer, # Required!
)

```

You can configure interrupt for agents and subagents on tool call as well as from within tool calls.
For more information, see [Human-in-the-loop](/oss/python/deepagents/human-in-the-loop).

## [​](#skills)Skills

You can use [skills](/oss/python/deepagents/overview) to provide your deep agent with new capabilities and expertise.
While [tools](/oss/python/deepagents/customization#tools) tend to cover lower level functionality like native file system actions or planning, skills can contain detailed instructions on how to complete tasks, reference info, and other assets, such as templates.
These files are only loaded by the agent when the agent has determined that the skill is useful for the current prompt.
This progressive disclosure reduces the amount of tokens and context the agent has to consider upon startup.
For example skills, see [Deep Agents example skills](https://github.com/langchain-ai/deepagentsjs/tree/main/examples/skills).
To add skills to your deep agent, pass them as an argument to `create_deep_agent`:

{' ' * (self.list_depth - 1)}- StateBackend
{' ' * (self.list_depth - 1)}- StoreBackend
{' ' * (self.list_depth - 1)}- FilesystemBackendGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from urllib.request import urlopen
from deepagents import create_deep_agent
from deepagents.backends import StateBackend
from deepagents.backends.utils import create_file_data
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.checkpoint.memory import MemorySaver

checkpointer = MemorySaver()
backend = StateBackend()

skill_url = "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/libs/cli/examples/skills/langgraph-docs/SKILL.md"
with urlopen(skill_url) as response:
 skill_content = response.read().decode('utf-8')

skills_files = {
 "/skills/langgraph-docs/SKILL.md": create_file_data(skill_content),
}

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 skills=["/skills/"],
 checkpointer=checkpointer,
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)], # for interpreter skills
)

result = agent.invoke(
 {
 "messages": [{"role": "user", "content": "What is langgraph?"}],
 # Seed the default StateBackend's in-state filesystem (virtual paths must start with "/").
 "files": skills_files,
 },
 config={"configurable": {"thread_id": "12345"}},
)

```

```
from urllib.request import urlopen
from deepagents import create_deep_agent
from deepagents.backends import StoreBackend
from deepagents.backends.utils import create_file_data
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.store.memory import InMemoryStore

store = InMemoryStore()
backend = StoreBackend(namespace=lambda _rt: ("filesystem",))

skill_url = "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/libs/cli/examples/skills/langgraph-docs/SKILL.md"
with urlopen(skill_url) as response:
 skill_content = response.read().decode('utf-8')

store.put(
 namespace=("filesystem",),
 key="/skills/langgraph-docs/SKILL.md",
 value=create_file_data(skill_content),
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 store=store,
 skills=["/skills/"],
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)],
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What is langgraph?"}]},
 config={"configurable": {"thread_id": "12345"}},
)

```

```
from deepagents import create_deep_agent
from deepagents.backends.filesystem import FilesystemBackend
from langchain_quickjs import CodeInterpreterMiddleware
from langgraph.checkpoint.memory import MemorySaver

# Checkpointer is REQUIRED for human-in-the-loop
checkpointer = MemorySaver()
root_dir = "/Users/user/{project}"
backend = FilesystemBackend(root_dir=root_dir)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 skills=[str(Path(root_dir) / "skills")],
 interrupt_on={
 "write_file": True,
 "read_file": False,
 "edit_file": True,
 },
 checkpointer=checkpointer, # Required!
 middleware=[CodeInterpreterMiddleware(skills_backend=backend)], # for interpreter skills
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What is langgraph?"}]},
 config={"configurable": {"thread_id": "12345"}},
)

```

## [​](#memory)Memory

Use [`AGENTS.md` files](https://agents.md/) to provide extra context to your deep agent.
You can pass one or more file paths to the `memory` parameter when creating your deep agent:

{' ' * (self.list_depth - 1)}- StateBackend
{' ' * (self.list_depth - 1)}- StoreBackend
{' ' * (self.list_depth - 1)}- FilesystemBackendGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from urllib.request import urlopen

from deepagents import create_deep_agent
from deepagents.backends.utils import create_file_data
from langgraph.checkpoint.memory import MemorySaver

with urlopen(
 "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/examples/text-to-sql-agent/AGENTS.md"
) as response:
 agents_md = response.read().decode("utf-8")
checkpointer = MemorySaver()

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 memory=[
 "/AGENTS.md"
 ],
 checkpointer=checkpointer,
)

result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Please tell me what's in your memory files.",
 }
 ],
 # Seed the default StateBackend's in-state filesystem (virtual paths must start with "/").
 "files": {"/AGENTS.md": create_file_data(agents_md)},
 },
 config={"configurable": {"thread_id": "123456"}},
)

```
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from urllib.request import urlopen

from deepagents import create_deep_agent
from deepagents.backends import StoreBackend
from deepagents.backends.utils import create_file_data
from langgraph.store.memory import InMemoryStore

with urlopen(
 "https://raw.githubusercontent.com/langchain-ai/deepagents/refs/heads/main/examples/text-to-sql-agent/AGENTS.md"
) as response:
 agents_md = response.read().decode("utf-8")

# Create the store and add the file to it
store = InMemoryStore()
file_data = create_file_data(agents_md)
store.put(
 namespace=("filesystem",),
 key="/AGENTS.md",
 value=file_data,
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=StoreBackend(namespace=lambda _rt: ("filesystem",)),
 store=store,
 memory=["/AGENTS.md"],
)

result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Please tell me what's in your memory files.",
 }
 ],
 "files": {"/AGENTS.md": create_file_data(agents_md)},
 },
 config={"configurable": {"thread_id": "12345"}},
)

```
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from deepagents import create_deep_agent
from deepagents.backends import FilesystemBackend
from langgraph.checkpoint.memory import MemorySaver

# Checkpointer is REQUIRED for human-in-the-loop
checkpointer = MemorySaver()

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=FilesystemBackend(root_dir="/Users/user/{project}"),
 memory=[
 "./AGENTS.md"
 ],
 interrupt_on={
 "write_file": True, # Default: approve, edit, reject
 "read_file": False, # No interrupts needed
 "edit_file": True, # Default: approve, edit, reject
 },
 checkpointer=checkpointer, # Required!
)

```

## [​](#profiles)Profiles

A [harness profile](/oss/python/deepagents/profiles#harness-profiles) is a reusable bundle of per-model configuration that `create_deep_agent` applies automatically when the matching model is selected. Profiles are the right tool when you want behaviour that follows the model—not the call site—such as a system prompt suffix tuned for Claude’s instruction style, tool descriptions rewritten for GPT, or extra middleware that only makes sense with a specific provider.
A single profile can carry: a custom base system prompt (`base_system_prompt`), an appended suffix (`system_prompt_suffix`), tool description overrides, tools or middleware to exclude, additional middleware to inject, and edits to the auto-added general-purpose subagent.

```
from deepagents import HarnessProfile, register_harness_profile

# Append a system-prompt suffix whenever gpt-5.4 is selected.
register_harness_profile(
 "openai:gpt-5.4",
 HarnessProfile(system_prompt_suffix="Respond in under 100 words."),
)

```

See [Profiles](/oss/python/deepagents/profiles) for registration keys, merge semantics, and plugin packaging. A narrower companion API, [provider profiles](/oss/python/deepagents/profiles#provider-profiles), packages model-construction arguments (API keys, timeouts, retry settings) for a provider.

## [​](#structured-output)Structured output

Deep Agents support [structured output](/oss/python/langchain/structured-output).
You can set a desired structured output schema by passing it as the `response_format` argument to the call to `create_deep_agent()`.
When the model generates the structured data, it’s captured, validated, and returned in the ‘structured_response’ key of the deep agent’s state.

```
import os
from typing import Literal

from pydantic import BaseModel, Field
from tavily import TavilyClient

from deepagents import create_deep_agent

tavily_client = TavilyClient(api_key=os.environ["TAVILY_API_KEY"])

def internet_search(
 query: str,
 max_results: int = 5,
 topic: Literal["general", "news", "finance"] = "general",
 include_raw_content: bool = False,
):
 """Run a web search"""
 return tavily_client.search(
 query,
 max_results=max_results,
 include_raw_content=include_raw_content,
 topic=topic,
 )

class WeatherReport(BaseModel):
 """A structured weather report with current conditions and forecast."""
 location: str = Field(description="The location for this weather report")
 temperature: float = Field(description="Current temperature in Celsius")
 condition: str = Field(
 description="Current weather condition (e.g., sunny, cloudy, rainy)"
 )
 humidity: int = Field(description="Humidity percentage")
 wind_speed: float = Field(description="Wind speed in km/h")
 forecast: str = Field(description="Brief forecast for the next 24 hours")

agent = create_deep_agent(
 model=model,
 response_format=WeatherReport,
 tools=[internet_search],
)

result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "What's the weather like in San Francisco?",
 }
 ]
 }
)

print(result["structured_response"])
# location='San Francisco, California' temperature=18.3 condition='Sunny' humidity=48 wind_speed=7.6 forecast='Pleasant sunny conditions expected to continue with temperatures around 64°F (18°C) during the day, dropping to around 52°F (11°C) at night. Clear skies with minimal precipitation expected.'

```

For more information and examples, see [response format](/oss/python/langchain/structured-output#response-format).

## [​](#advanced)Advanced

`create_deep_agent` pre-assembles a middleware stack on top of [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent). To build a fully custom agent—choosing exactly which capabilities to include—see [Configure the harness](/oss/python/langchain/agents#configure-the-harness).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/customization.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[QuickstartPrevious](/oss/python/deepagents/quickstart)[Comparison with Claude Agent SDKNext](/oss/python/deepagents/comparison)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
