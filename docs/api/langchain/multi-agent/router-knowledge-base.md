# Build a multi-source knowledge base with routing

[Tutorials](/oss/python/deepagents/data-analysis)[Multi-agent](/oss/python/langchain/multi-agent/subagents-personal-assistant)

# Build a multi-source knowledge base with routing
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

The **router pattern** is a [multi-agent](/oss/python/langchain/multi-agent) architecture where a routing step classifies input and directs it to specialized agents, with results synthesized into a combined response. This pattern excels when your organization’s knowledge lives across distinct **verticals** (separate knowledge domains that each require their own agent with specialized tools and prompts).
In this tutorial, you’ll build a multi-source knowledge base router that demonstrates these benefits through a realistic enterprise scenario. The system will coordinate three specialists:

{' ' * (self.list_depth - 1)}- A **GitHub agent** that searches code, issues, and pull requests.

{' ' * (self.list_depth - 1)}- A **Notion agent** that searches internal documentation and wikis.

{' ' * (self.list_depth - 1)}- A **Slack agent** that searches relevant threads and discussions.

When a user asks “How do I authenticate API requests?”, the router decomposes the query into source-specific sub-questions, routes them to the relevant agents in parallel, and synthesizes results into a coherent answer.

### [​](#why-use-a-router)Why use a router?

The router pattern provides several advantages:

{' ' * (self.list_depth - 1)}- **Parallel execution**: Query multiple sources simultaneously, reducing latency compared to sequential approaches.

{' ' * (self.list_depth - 1)}- **Specialized agents**: Each vertical has focused tools and prompts optimized for its domain.

{' ' * (self.list_depth - 1)}- **Selective routing**: Not every query needs every source—the router intelligently selects relevant verticals.

{' ' * (self.list_depth - 1)}- **Targeted sub-questions**: Each agent receives a question tailored to its domain, improving result quality.

{' ' * (self.list_depth - 1)}- **Clean synthesis**: Results from multiple sources are combined into a single, coherent response.

### [​](#concepts)Concepts

We will cover the following concepts:

{' ' * (self.list_depth - 1)}- [Multi-agent systems](/oss/python/langchain/multi-agent)

{' ' * (self.list_depth - 1)}- [StateGraph](/oss/python/langgraph/graph-api) for workflow orchestration

{' ' * (self.list_depth - 1)}- [Send API](/oss/python/langgraph/graph-api#send) for parallel execution

**Router vs. Subagents**: The [subagents pattern](/oss/python/langchain/multi-agent/subagents) can also route to multiple agents. Use the router pattern when you need specialized preprocessing, custom routing logic, or want explicit control over parallel execution. Use the subagents pattern when you want the LLM to decide which agents to call dynamically.

## [​](#setup)Setup

### [​](#installation)Installation

This tutorial requires the `langchain` and `langgraph` packages:
pipuvconda
```
pip install langchain langgraph

```

For more details, see our [Installation guide](/oss/python/langchain/install).

### [​](#langsmith)LangSmith

Set up [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-router-knowledge-base) to inspect what is happening inside your agent. Then set the following environment variables:
bashpython
```
export LANGSMITH_TRACING="true"
export LANGSMITH_API_KEY="..."

```

### [​](#select-an-llm)Select an LLM

Select a chat model from LangChain’s suite of integrations:

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- Azure
{' ' * (self.list_depth - 1)}- Google Gemini
{' ' * (self.list_depth - 1)}- AWS Bedrock
{' ' * (self.list_depth - 1)}- HuggingFace
{' ' * (self.list_depth - 1)}- OpenRouter👉 Read the [OpenAI chat model integration docs](/oss/python/integrations/chat/openai)
```
pip install -U "langchain[openai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["OPENAI_API_KEY"] = "sk-..."

model = init_chat_model("gpt-5.4")

```
👉 Read the [Anthropic chat model integration docs](/oss/python/integrations/chat/anthropic)
```
pip install -U "langchain[anthropic]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["ANTHROPIC_API_KEY"] = "sk-..."

model = init_chat_model("claude-sonnet-4-6")

```
👉 Read the [Azure chat model integration docs](/oss/python/integrations/chat/azure_chat_openai)
```
pip install -U "langchain[openai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["AZURE_OPENAI_API_KEY"] = "..."
os.environ["AZURE_OPENAI_ENDPOINT"] = "..."
os.environ["OPENAI_API_VERSION"] = "2025-03-01-preview"

model = init_chat_model(
 "azure_openai:gpt-5.4",
 azure_deployment=os.environ["AZURE_OPENAI_DEPLOYMENT_NAME"],
)

```
👉 Read the [Google GenAI chat model integration docs](/oss/python/integrations/chat/google_generative_ai)
```
pip install -U "langchain[google-genai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["GOOGLE_API_KEY"] = "..."

model = init_chat_model("google_genai:gemini-2.5-flash-lite")

```
👉 Read the [AWS Bedrock chat model integration docs](/oss/python/integrations/chat/bedrock)
```
pip install -U "langchain[aws]"

```
init_chat_modelModel Class
```
from langchain.chat_models import init_chat_model

# Follow the steps here to configure your credentials:
# https://docs.aws.amazon.com/bedrock/latest/userguide/getting-started.html

model = init_chat_model(
 "anthropic.claude-3-5-sonnet-20240620-v1:0",
 model_provider="bedrock_converse",
)

```
👉 Read the [HuggingFace chat model integration docs](/oss/python/integrations/chat/huggingface)
```
pip install -U "langchain[huggingface]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["HUGGINGFACEHUB_API_TOKEN"] = "hf_..."

model = init_chat_model(
 "microsoft/Phi-3-mini-4k-instruct",
 model_provider="huggingface",
 temperature=0.7,
 max_tokens=1024,
)

```
👉 Read the [OpenRouter chat model integration docs](/oss/python/integrations/chat/openrouter)
```
pip install -U "langchain-openrouter"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["OPENROUTER_API_KEY"] = "sk-..."

model = init_chat_model(
 "auto",
 model_provider="openrouter",
)

```

## [​](#1-define-state)1. Define state

First, define the state schemas. We use three types:

{' ' * (self.list_depth - 1)}- **`AgentInput`**: Simple state passed to each subagent (just a query)

{' ' * (self.list_depth - 1)}- **`AgentOutput`**: Result returned by each subagent (source name + result)

{' ' * (self.list_depth - 1)}- **`RouterState`**: Main workflow state tracking the query, classifications, results, and final answer

```
from typing import Annotated, Literal, TypedDict
import operator

class AgentInput(TypedDict):
 """Simple input state for each subagent."""
 query: str

class AgentOutput(TypedDict):
 """Output from each subagent."""
 source: str
 result: str

class Classification(TypedDict):
 """A single routing decision: which agent to call with what query."""
 source: Literal["github", "notion", "slack"]
 query: str

class RouterState(TypedDict):
 query: str
 classifications: list[Classification]
 results: Annotated[list[AgentOutput], operator.add] # Reducer collects parallel results
 final_answer: str

```

The `results` field uses a **reducer** (`operator.add` in Python, a concat function in JS) to collect outputs from parallel agent executions into a single list.

## [​](#2-define-tools-for-each-vertical)2. Define tools for each vertical

Create tools for each knowledge domain. In a production system, these would call actual APIs. For this tutorial, we use stub implementations that return mock data. We define 7 tools across 3 verticals: GitHub (search code, issues, PRs), Notion (search docs, get page), and Slack (search messages, get thread).

```
from langchain.tools import tool

@tool
def search_code(query: str, repo: str = "main") -> str:
 """Search code in GitHub repositories."""
 return f"Found code matching '{query}' in {repo}: authentication middleware in src/auth.py"

@tool
def search_issues(query: str) -> str:
 """Search GitHub issues and pull requests."""
 return f"Found 3 issues matching '{query}': #142 (API auth docs), #89 (OAuth flow), #203 (token refresh)"

@tool
def search_prs(query: str) -> str:
 """Search pull requests for implementation details."""
 return f"PR #156 added JWT authentication, PR #178 updated OAuth scopes"

@tool
def search_notion(query: str) -> str:
 """Search Notion workspace for documentation."""
 return f"Found documentation: 'API Authentication Guide' - covers OAuth2 flow, API keys, and JWT tokens"

@tool
def get_page(page_id: str) -> str:
 """Get a specific Notion page by ID."""
 return f"Page content: Step-by-step authentication setup instructions"

@tool
def search_slack(query: str) -> str:
 """Search Slack messages and threads."""
 return f"Found discussion in #engineering: 'Use Bearer tokens for API auth, see docs for refresh flow'"

@tool
def get_thread(thread_id: str) -> str:
 """Get a specific Slack thread."""
 return f"Thread discusses best practices for API key rotation"

```
See all 43 lines

## [​](#3-create-specialized-agents)3. Create specialized agents

Create an agent for each vertical. Each agent has domain-specific tools and a prompt optimized for its knowledge source. All three follow the same pattern—only the tools and system prompt differ.

```
from langchain.agents import create_agent
from langchain.chat_models import init_chat_model

model = init_chat_model("openai:gpt-5.4")

github_agent = create_agent(
 model,
 tools=[search_code, search_issues, search_prs],
 system_prompt=(
 "You are a GitHub expert. Answer questions about code, "
 "API references, and implementation details by searching "
 "repositories, issues, and pull requests."
 ),
)

notion_agent = create_agent(
 model,
 tools=[search_notion, get_page],
 system_prompt=(
 "You are a Notion expert. Answer questions about internal "
 "processes, policies, and team documentation by searching "
 "the organization's Notion workspace."
 ),
)

slack_agent = create_agent(
 model,
 tools=[search_slack, get_thread],
 system_prompt=(
 "You are a Slack expert. Answer questions by searching "
 "relevant threads and discussions where team members have "
 "shared knowledge and solutions."
 ),
)

```
See all 34 lines

## [​](#4-build-the-router-workflow)4. Build the router workflow

Now build the router workflow using a StateGraph. The workflow has four main steps:

{' ' * (self.list_depth - 1)}- **Classify**: Analyze the query and determine which agents to invoke with what sub-questions

{' ' * (self.list_depth - 1)}- **Route**: Fan out to selected agents in parallel using `Send`

{' ' * (self.list_depth - 1)}- **Query agents**: Each agent receives a simple `AgentInput` and returns an `AgentOutput`

{' ' * (self.list_depth - 1)}- **Synthesize**: Combine collected results into a coherent response

```
from pydantic import BaseModel, Field
from langgraph.graph import StateGraph, START, END
from langgraph.types import Send

router_llm = init_chat_model("openai:gpt-5.4-mini")

# Define structured output schema for the classifier
class ClassificationResult(BaseModel):
 """Result of classifying a user query into agent-specific sub-questions."""
 classifications: list[Classification] = Field(
 description="List of agents to invoke with their targeted sub-questions"
 )

def classify_query(state: RouterState) -> dict:
 """Classify query and determine which agents to invoke."""
 structured_llm = router_llm.with_structured_output(ClassificationResult)

 result = structured_llm.invoke([
 {
 "role": "system",
 "content": """Analyze this query and determine which knowledge bases to consult.
For each relevant source, generate a targeted sub-question optimized for that source.

Available sources:
- github: Code, API references, implementation details, issues, pull requests
- notion: Internal documentation, processes, policies, team wikis
- slack: Team discussions, informal knowledge sharing, recent conversations

Return ONLY the sources that are relevant to the query. Each source should have
a targeted sub-question optimized for that specific knowledge domain.

Example for "How do I authenticate API requests?":
- github: "What authentication code exists? Search for auth middleware, JWT handling"
- notion: "What authentication documentation exists? Look for API auth guides"
(slack omitted because it's not relevant for this technical question)"""
 },
 {"role": "user", "content": state["query"]}
 ])

 return {"classifications": result.classifications}

def route_to_agents(state: RouterState) -> list[Send]:
 """Fan out to agents based on classifications."""
 return [
 Send(c["source"], {"query": c["query"]})
 for c in state["classifications"]
 ]

def query_github(state: AgentInput) -> dict:
 """Query the GitHub agent."""
 result = github_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "github", "result": result["messages"][-1].content}]}

def query_notion(state: AgentInput) -> dict:
 """Query the Notion agent."""
 result = notion_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "notion", "result": result["messages"][-1].content}]}

def query_slack(state: AgentInput) -> dict:
 """Query the Slack agent."""
 result = slack_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "slack", "result": result["messages"][-1].content}]}

def synthesize_results(state: RouterState) -> dict:
 """Combine results from all agents into a coherent answer."""
 if not state["results"]:
 return {"final_answer": "No results found from any knowledge source."}

 # Format results for synthesis
 formatted = [
 f"**From {r['source'].title()}:**\n{r['result']}"
 for r in state["results"]
 ]

 synthesis_response = router_llm.invoke([
 {
 "role": "system",
 "content": f"""Synthesize these search results to answer the original question: "{state['query']}"

- Combine information from multiple sources without redundancy
- Highlight the most relevant and actionable information
- Note any discrepancies between sources
- Keep the response concise and well-organized"""
 },
 {"role": "user", "content": "\n\n".join(formatted)}
 ])

 return {"final_answer": synthesis_response.content}

```

## [​](#5-compile-the-workflow)5. Compile the workflow

Now assemble the workflow by connecting nodes with edges. The key is using `add_conditional_edges` with the routing function to enable parallel execution:

```
workflow = (
 StateGraph(RouterState)
 .add_node("classify", classify_query)
 .add_node("github", query_github)
 .add_node("notion", query_notion)
 .add_node("slack", query_slack)
 .add_node("synthesize", synthesize_results)
 .add_edge(START, "classify")
 .add_conditional_edges("classify", route_to_agents, ["github", "notion", "slack"])
 .add_edge("github", "synthesize")
 .add_edge("notion", "synthesize")
 .add_edge("slack", "synthesize")
 .add_edge("synthesize", END)
 .compile()
)

```

The `add_conditional_edges` call connects the classify node to the agent nodes through the `route_to_agents` function. When `route_to_agents` returns multiple `Send` objects, those nodes execute in parallel.

## [​](#6-use-the-router)6. Use the router

Test your router with queries that span multiple knowledge domains:

```
result = workflow.invoke({
 "query": "How do I authenticate API requests?"
})

print("Original query:", result["query"])
print("\nClassifications:")
for c in result["classifications"]:
 print(f" {c['source']}: {c['query']}")
print("\n" + "=" * 60 + "\n")
print("Final Answer:")
print(result["final_answer"])

```

Expected output:

```
Original query: How do I authenticate API requests?

Classifications:
 github: What authentication code exists? Search for auth middleware, JWT handling
 notion: What authentication documentation exists? Look for API auth guides

============================================================

Final Answer:
To authenticate API requests, you have several options:

1. **JWT Tokens**: The recommended approach for most use cases.
 Implementation details are in `src/auth.py` (PR #156).

2. **OAuth2 Flow**: For third-party integrations, follow the OAuth2
 flow documented in Notion's 'API Authentication Guide'.

3. **API Keys**: For server-to-server communication, use Bearer tokens
 in the Authorization header.

For token refresh handling, see issue #203 and PR #178 for the latest
OAuth scope updates.

```

The router analyzed the query, classified it to determine which agents to invoke (GitHub and Notion, but not Slack for this technical question), queried both agents in parallel, and synthesized the results into a coherent answer.

## [​](#7-understanding-the-architecture)7. Understanding the architecture

The router workflow follows a clear pattern:

### [​](#classification-phase)Classification phase

The `classify_query` function uses **structured output** to analyze the user’s query and determine which agents to invoke. This is where the routing intelligence lives:

{' ' * (self.list_depth - 1)}- Uses a Pydantic model (Python) or Zod schema (JS) to ensure valid output

{' ' * (self.list_depth - 1)}- Returns a list of `Classification` objects, each with a `source` and targeted `query`

{' ' * (self.list_depth - 1)}- Only includes relevant sources—irrelevant ones are simply omitted

This structured approach is more reliable than free-form JSON parsing and makes the routing logic explicit.

### [​](#parallel-execution-with-send)Parallel execution with send

The `route_to_agents` function maps classifications to `Send` objects. Each `Send` specifies the target node and the state to pass:

```
# Classifications: [{"source": "github", "query": "..."}, {"source": "notion", "query": "..."}]
# Becomes:
[Send("github", {"query": "..."}), Send("notion", {"query": "..."})]
# Both agents execute simultaneously, each receiving only the query it needs

```

Each agent node receives a simple `AgentInput` with just a `query` field—not the full router state. This keeps the interface clean and explicit.

### [​](#result-collection-with-reducers)Result collection with reducers

Agent results flow back to the main state via a **reducer**. Each agent returns:

```
{"results": [{"source": "github", "result": "..."}]}

```

The reducer (`operator.add` in Python) concatenates these lists, collecting all parallel results into `state["results"]`.

### [​](#synthesis-phase)Synthesis phase

After all agents complete, the `synthesize_results` function iterates over the collected results:

{' ' * (self.list_depth - 1)}- Waits for all parallel branches to complete (LangGraph handles this automatically)

{' ' * (self.list_depth - 1)}- References the original query to ensure the answer addresses what the user asked

{' ' * (self.list_depth - 1)}- Combines information from all sources without redundancy

**Partial results**: In this tutorial, all selected agents must complete before synthesis.

## [​](#8-complete-working-example)8. Complete working example

Here’s everything together in a runnable script:

Show View complete code
```
"""
Multi-Source Knowledge Router Example

This example demonstrates the router pattern for multi-agent systems.
A router classifies queries, routes them to specialized agents in parallel,
and synthesizes results into a combined response.
"""

import operator
from typing import Annotated, Literal, TypedDict

from langchain.agents import create_agent
from langchain.chat_models import init_chat_model
from langchain.tools import tool
from langgraph.graph import StateGraph, START, END
from langgraph.types import Send
from pydantic import BaseModel, Field

# State definitions
class AgentInput(TypedDict):
 """Simple input state for each subagent."""
 query: str

class AgentOutput(TypedDict):
 """Output from each subagent."""
 source: str
 result: str

class Classification(TypedDict):
 """A single routing decision: which agent to call with what query."""
 source: Literal["github", "notion", "slack"]
 query: str

class RouterState(TypedDict):
 query: str
 classifications: list[Classification]
 results: Annotated[list[AgentOutput], operator.add]
 final_answer: str

# Structured output schema for classifier
class ClassificationResult(BaseModel):
 """Result of classifying a user query into agent-specific sub-questions."""
 classifications: list[Classification] = Field(
 description="List of agents to invoke with their targeted sub-questions"
 )

# Tools
@tool
def search_code(query: str, repo: str = "main") -> str:
 """Search code in GitHub repositories."""
 return f"Found code matching '{query}' in {repo}: authentication middleware in src/auth.py"

@tool
def search_issues(query: str) -> str:
 """Search GitHub issues and pull requests."""
 return f"Found 3 issues matching '{query}': #142 (API auth docs), #89 (OAuth flow), #203 (token refresh)"

@tool
def search_prs(query: str) -> str:
 """Search pull requests for implementation details."""
 return f"PR #156 added JWT authentication, PR #178 updated OAuth scopes"

@tool
def search_notion(query: str) -> str:
 """Search Notion workspace for documentation."""
 return f"Found documentation: 'API Authentication Guide' - covers OAuth2 flow, API keys, and JWT tokens"

@tool
def get_page(page_id: str) -> str:
 """Get a specific Notion page by ID."""
 return f"Page content: Step-by-step authentication setup instructions"

@tool
def search_slack(query: str) -> str:
 """Search Slack messages and threads."""
 return f"Found discussion in #engineering: 'Use Bearer tokens for API auth, see docs for refresh flow'"

@tool
def get_thread(thread_id: str) -> str:
 """Get a specific Slack thread."""
 return f"Thread discusses best practices for API key rotation"

# Models and agents
model = init_chat_model("openai:gpt-5.4")
router_llm = init_chat_model("openai:gpt-5.4-mini")

github_agent = create_agent(
 model,
 tools=[search_code, search_issues, search_prs],
 system_prompt=(
 "You are a GitHub expert. Answer questions about code, "
 "API references, and implementation details by searching "
 "repositories, issues, and pull requests."
 ),
)

notion_agent = create_agent(
 model,
 tools=[search_notion, get_page],
 system_prompt=(
 "You are a Notion expert. Answer questions about internal "
 "processes, policies, and team documentation by searching "
 "the organization's Notion workspace."
 ),
)

slack_agent = create_agent(
 model,
 tools=[search_slack, get_thread],
 system_prompt=(
 "You are a Slack expert. Answer questions by searching "
 "relevant threads and discussions where team members have "
 "shared knowledge and solutions."
 ),
)

# Workflow nodes
def classify_query(state: RouterState) -> dict:
 """Classify query and determine which agents to invoke."""
 structured_llm = router_llm.with_structured_output(ClassificationResult)

 result = structured_llm.invoke([
 {
 "role": "system",
 "content": """Analyze this query and determine which knowledge bases to consult.
For each relevant source, generate a targeted sub-question optimized for that source.

Available sources:
- github: Code, API references, implementation details, issues, pull requests
- notion: Internal documentation, processes, policies, team wikis
- slack: Team discussions, informal knowledge sharing, recent conversations

Return ONLY the sources that are relevant to the query."""
 },
 {"role": "user", "content": state["query"]}
 ])

 return {"classifications": result.classifications}

def route_to_agents(state: RouterState) -> list[Send]:
 """Fan out to agents based on classifications."""
 return [
 Send(c["source"], {"query": c["query"]})
 for c in state["classifications"]
 ]

def query_github(state: AgentInput) -> dict:
 """Query the GitHub agent."""
 result = github_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "github", "result": result["messages"][-1].content}]}

def query_notion(state: AgentInput) -> dict:
 """Query the Notion agent."""
 result = notion_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "notion", "result": result["messages"][-1].content}]}

def query_slack(state: AgentInput) -> dict:
 """Query the Slack agent."""
 result = slack_agent.invoke({
 "messages": [{"role": "user", "content": state["query"]}]
 })
 return {"results": [{"source": "slack", "result": result["messages"][-1].content}]}

def synthesize_results(state: RouterState) -> dict:
 """Combine results from all agents into a coherent answer."""
 if not state["results"]:
 return {"final_answer": "No results found from any knowledge source."}

 formatted = [
 f"**From {r['source'].title()}:**\n{r['result']}"
 for r in state["results"]
 ]

 synthesis_response = router_llm.invoke([
 {
 "role": "system",
 "content": f"""Synthesize these search results to answer the original question: "{state['query']}"

- Combine information from multiple sources without redundancy
- Highlight the most relevant and actionable information
- Note any discrepancies between sources
- Keep the response concise and well-organized"""
 },
 {"role": "user", "content": "\n\n".join(formatted)}
 ])

 return {"final_answer": synthesis_response.content}

# Build workflow
workflow = (
 StateGraph(RouterState)
 .add_node("classify", classify_query)
 .add_node("github", query_github)
 .add_node("notion", query_notion)
 .add_node("slack", query_slack)
 .add_node("synthesize", synthesize_results)
 .add_edge(START, "classify")
 .add_conditional_edges("classify", route_to_agents, ["github", "notion", "slack"])
 .add_edge("github", "synthesize")
 .add_edge("notion", "synthesize")
 .add_edge("slack", "synthesize")
 .add_edge("synthesize", END)
 .compile()
)

if __name__ == "__main__":
 result = workflow.invoke({
 "query": "How do I authenticate API requests?"
 })

 print("Original query:", result["query"])
 print("\nClassifications:")
 for c in result["classifications"]:
 print(f" {c['source']}: {c['query']}")
 print("\n" + "=" * 60 + "\n")
 print("Final Answer:")
 print(result["final_answer"])

```

## [​](#9-advanced-stateful-routers)9. Advanced: Stateful routers

The router we’ve built so far is **stateless** (each request is handled independently with no memory between calls). For multi-turn conversations, you need a **stateful** approach.

### [​](#tool-wrapper-approach)Tool wrapper approach

The simplest way to add conversation memory is to wrap the stateless router as a tool that a conversational agent can call:

```
from langgraph.checkpoint.memory import InMemorySaver

@tool
def search_knowledge_base(query: str) -> str:
 """Search across multiple knowledge sources (GitHub, Notion, Slack).

 Use this to find information about code, documentation, or team discussions.
 """
 result = workflow.invoke({"query": query})
 return result["final_answer"]

conversational_agent = create_agent(
 model,
 tools=[search_knowledge_base],
 system_prompt=(
 "You are a helpful assistant that answers questions about our organization. "
 "Use the search_knowledge_base tool to find information across our code, "
 "documentation, and team discussions."
 ),
 checkpointer=InMemorySaver(),
)

```

This approach keeps the router stateless while the conversational agent handles memory and context. The user can have a multi-turn conversation, and the agent will call the router tool as needed.

```
config = {"configurable": {"thread_id": "user-123"}}

result = conversational_agent.invoke(
 {"messages": [{"role": "user", "content": "How do I authenticate API requests?"}]},
 config
)
print(result["messages"][-1].content)

result = conversational_agent.invoke(
 {"messages": [{"role": "user", "content": "What about rate limiting for those endpoints?"}]},
 config
)
print(result["messages"][-1].content)

```

The tool wrapper approach is recommended for most use cases. It provides clean separation: the router handles multi-source querying, while the conversational agent handles context and memory.

### [​](#full-persistence-approach)Full persistence approach

If you need the router itself to maintain state—for example, to use previous search results in routing decisions—use [persistence](/oss/python/langchain/short-term-memory) to store message history at the router level.
**Stateful routers add complexity.** When routing to different agents across turns, conversations may feel inconsistent if agents have different tones or prompts. Consider the [handoffs pattern](/oss/python/langchain/multi-agent/handoffs) or [subagents pattern](/oss/python/langchain/multi-agent/subagents) instead—both provide clearer semantics for multi-turn conversations with different agents.

## [​](#10-key-takeaways)10. Key takeaways

The router pattern excels when you have:

{' ' * (self.list_depth - 1)}- **Distinct verticals**: Separate knowledge domains that each require specialized tools and prompts

{' ' * (self.list_depth - 1)}- **Parallel query needs**: Questions that benefit from querying multiple sources simultaneously

{' ' * (self.list_depth - 1)}- **Synthesis requirements**: Results from multiple sources need to be combined into a coherent response

The pattern has three phases: **decompose** (analyze the query and generate targeted sub-questions), **route** (execute queries in parallel), and **synthesize** (combine results).
**When to use the router pattern**Use the router pattern when you have multiple independent knowledge sources, need low-latency parallel queries, and want explicit control over routing logic.For simpler cases with dynamic tool selection, consider the [subagents pattern](/oss/python/langchain/multi-agent/subagents). For workflows where agents need to converse with users sequentially, consider [handoffs](/oss/python/langchain/multi-agent/handoffs).

## [​](#next-steps)Next steps

{' ' * (self.list_depth - 1)}- Learn about [handoffs](/oss/python/langchain/multi-agent/handoffs) for agent-to-agent conversations

{' ' * (self.list_depth - 1)}- Explore the [subagents pattern](/oss/python/langchain/multi-agent/subagents-personal-assistant) for centralized orchestration

{' ' * (self.list_depth - 1)}- Read the [multi-agent overview](/oss/python/langchain/multi-agent) to compare different patterns

{' ' * (self.list_depth - 1)}- Use [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-router-knowledge-base) to debug and monitor your router

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/multi-agent/router-knowledge-base.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build customer support with handoffsPrevious](/oss/python/langchain/multi-agent/handoffs-customer-support)[Build a SQL assistant with on-demand skillsNext](/oss/python/langchain/multi-agent/skills-sql-assistant)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
