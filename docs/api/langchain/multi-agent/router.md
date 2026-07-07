# Router

[Advanced usage](/oss/python/langchain/guardrails)[Multi-agent](/oss/python/langchain/multi-agent/index)

# Router
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.In the **router** architecture, a routing step classifies input and directs it to specialized [agents](/oss/python/langchain/agents). This is useful when you have distinct **verticals** (separate knowledge domains that each require their own agent).

## [​](#key-characteristics)Key characteristics

{' ' * (self.list_depth - 1)}- Router decomposes the query

{' ' * (self.list_depth - 1)}- Zero or more specialized agents are invoked in parallel

{' ' * (self.list_depth - 1)}- Results are synthesized into a coherent response

## [​](#when-to-use)When to use

Use the router pattern when you have distinct verticals (separate knowledge domains that each require their own agent), need to query multiple sources in parallel, and want to synthesize results into a combined response.

## [​](#basic-implementation)Basic implementation

The router classifies the query and directs it to the appropriate agent(s). Use [`Command`](/oss/python/langgraph/graph-api#command) for single-agent routing or [`Send`](/oss/python/langgraph/graph-api#send) for parallel fan-out to multiple agents.

{' ' * (self.list_depth - 1)}- Single agent
{' ' * (self.list_depth - 1)}- Multiple agents (parallel)Use `Command` to route to a single specialized agent:
```
from langgraph.types import Command

def classify_query(query: str) -> str:
 """Use LLM to classify query and determine the appropriate agent."""
 # Classification logic here
 ...

def route_query(state: State) -> Command:
 """Route to the appropriate agent based on query classification."""
 active_agent = classify_query(state["query"])

 # Route to the selected agent
 return Command(goto=active_agent)

```
Use `Send` to fan out to multiple specialized agents in parallel:
```
from typing import TypedDict
from langgraph.types import Send

class ClassificationResult(TypedDict):
 query: str
 agent: str

def classify_query(query: str) -> list[ClassificationResult]:
 """Use LLM to classify query and determine which agents to invoke."""
 # Classification logic here
 ...

def route_query(state: State):
 """Route to relevant agents based on query classification."""
 classifications = classify_query(state["query"])

 # Fan out to selected agents in parallel
 return [
 Send(c["agent"], {"query": c["query"]})
 for c in classifications
 ]

```

For a complete implementation, see the tutorial below.

## Tutorial: Build a multi-source knowledge base with routing
Build a router that queries GitHub, Notion, and Slack in parallel, then synthesizes results into a coherent answer. Covers state definition, specialized agents, parallel execution with `Send`, and result synthesis.

## [​](#stateless-vs-stateful)Stateless vs. stateful

Two approaches:

{' ' * (self.list_depth - 1)}- [**Stateless routers**](#stateless) address each request independently

{' ' * (self.list_depth - 1)}- [**Stateful routers**](#stateful) maintain conversation history across requests

## [​](#stateless)Stateless

Each request is routed independently—no memory between calls. For multi-turn conversations, see [Stateful routers](#stateful).
**Router vs. Subagents**: Both patterns can dispatch work to multiple agents, but they differ in how routing decisions are made:

{' ' * (self.list_depth - 1)}- **Router**: A dedicated routing step (often a single LLM call or rule-based logic) that classifies the input and dispatches to agents. The router itself typically doesn’t maintain conversation history or perform multi-turn orchestration—it’s a preprocessing step.

{' ' * (self.list_depth - 1)}- **Subagents**: An main supervisor agent dynamically decides which [subagents](/oss/python/langchain/multi-agent/subagents) to call as part of an ongoing conversation. The main agent maintains context, can call multiple subagents across turns, and orchestrates complex multi-step workflows.
Use a **router** when you have clear input categories and want deterministic or lightweight classification. Use a **supervisor** when you need flexible, conversation-aware orchestration where the LLM decides what to do next based on evolving context.

## [​](#stateful)Stateful

For multi-turn conversations, you need to maintain context across invocations.

### [​](#tool-wrapper)Tool wrapper

The simplest approach: wrap the stateless router as a tool that a conversational agent can call. The conversational agent handles memory and context; the router stays stateless. This avoids the complexity of managing conversation history across multiple parallel agents.

```
@tool
def search_docs(query: str) -> str:
 """Search across multiple documentation sources."""
 result = workflow.invoke({"query": query})
 return result["final_answer"]

# Conversational agent uses the router as a tool
conversational_agent = create_agent(
 model,
 tools=[search_docs],
 prompt="You are a helpful assistant. Use search_docs to answer questions."
)

```

### [​](#full-persistence)Full persistence

If you need the router itself to maintain state, use [persistence](/oss/python/langchain/short-term-memory) to store message history. When routing to an agent, fetch previous messages from state and selectively include them in the agent’s context—this is a lever for [context engineering](/oss/python/langchain/context-engineering).
**Stateful routers require custom history management.** If the router switches between agents across turns, conversations may not feel fluid to end users when agents have different tones or prompts. With parallel invocation, you’ll need to maintain history at the router level (inputs and synthesized outputs) and leverage this history in routing logic. Consider the [handoffs pattern](/oss/python/langchain/multi-agent/handoffs) or [subagents pattern](/oss/python/langchain/multi-agent/subagents) instead—both provide clearer semantics for multi-turn conversations.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/multi-agent/router.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[SkillsPrevious](/oss/python/langchain/multi-agent/skills)[Custom workflowNext](/oss/python/langchain/multi-agent/custom-workflow)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
