# Test

[Production](/oss/python/langgraph/application-structure)

# Test
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.After you’ve prototyped your LangGraph agent, a natural next step is to add tests. This guide covers some useful patterns you can use when writing unit tests.
Note that this guide is LangGraph-specific and covers scenarios around graphs with custom structures - if you are just getting started, check out [Test](/oss/python/langchain/test) that uses LangChain’s built-in [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) instead.

## [​](#prerequisites)Prerequisites

First, make sure you have [`pytest`](https://docs.pytest.org/) installed:

```
$ pip install -U pytest

```

## [​](#getting-started)Getting started

Because many LangGraph agents depend on state, a useful pattern is to create your graph before each test where you use it, then compile it within tests with a new checkpointer instance.
The below example shows how this works with a simple, linear graph that progresses through `node1` and `node2`. Each node updates the single state key `my_key`:

```
import pytest

from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

def create_graph() -> StateGraph:
 class MyState(TypedDict):
 my_key: str

 graph = StateGraph(MyState)
 graph.add_node("node1", lambda state: {"my_key": "hello from node1"})
 graph.add_node("node2", lambda state: {"my_key": "hello from node2"})
 graph.add_edge(START, "node1")
 graph.add_edge("node1", "node2")
 graph.add_edge("node2", END)
 return graph

def test_basic_agent_execution() -> None:
 checkpointer = MemorySaver()
 graph = create_graph()
 compiled_graph = graph.compile(checkpointer=checkpointer)
 result = compiled_graph.invoke(
 {"my_key": "initial_value"},
 config={"configurable": {"thread_id": "1"}}
 )
 assert result["my_key"] == "hello from node2"

```

## [​](#testing-individual-nodes-and-edges)Testing individual nodes and edges

Compiled LangGraph agents expose references to each individual node as `graph.nodes`. You can take advantage of this to test individual nodes within your agent. Note that this will bypass any checkpointers passed when compiling the graph:

```
import pytest

from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

def create_graph() -> StateGraph:
 class MyState(TypedDict):
 my_key: str

 graph = StateGraph(MyState)
 graph.add_node("node1", lambda state: {"my_key": "hello from node1"})
 graph.add_node("node2", lambda state: {"my_key": "hello from node2"})
 graph.add_edge(START, "node1")
 graph.add_edge("node1", "node2")
 graph.add_edge("node2", END)
 return graph

def test_individual_node_execution() -> None:
 # Will be ignored in this example
 checkpointer = MemorySaver()
 graph = create_graph()
 compiled_graph = graph.compile(checkpointer=checkpointer)
 # Only invoke node 1
 result = compiled_graph.nodes["node1"].invoke(
 {"my_key": "initial_value"},
 )
 assert result["my_key"] == "hello from node1"

```

## [​](#partial-execution)Partial execution

For agents made up of larger graphs, you may wish to test partial execution paths within your agent rather than the entire flow end-to-end. In some cases, it may make semantic sense to [restructure these sections as subgraphs](/oss/python/langgraph/use-subgraphs), which you can invoke in isolation as normal.
However, if you do not wish to make changes to your agent graph’s overall structure, you can use LangGraph’s persistence mechanisms to simulate a state where your agent is paused right before the beginning of the desired section, and will pause again at the end of the desired section. The steps are as follows:

{' ' * (self.list_depth - 1)}- Compile your agent with a checkpointer (the in-memory checkpointer [`InMemorySaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.memory.InMemorySaver) will suffice for testing).

{' ' * (self.list_depth - 1)}- Call your agent’s [`update_state`](/oss/python/langgraph/use-time-travel) method with an [`as_node`](/oss/python/langgraph/use-time-travel#from-a-specific-node) parameter set to the name of the node *before* the one you want to start your test.

{' ' * (self.list_depth - 1)}- Invoke your agent with the same `thread_id` you used to update the state and an `interrupt_after` parameter set to the name of the node you want to stop at.

Here’s an example that executes only the second and third nodes in a linear graph:

```
import pytest

from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

def create_graph() -> StateGraph:
 class MyState(TypedDict):
 my_key: str

 graph = StateGraph(MyState)
 graph.add_node("node1", lambda state: {"my_key": "hello from node1"})
 graph.add_node("node2", lambda state: {"my_key": "hello from node2"})
 graph.add_node("node3", lambda state: {"my_key": "hello from node3"})
 graph.add_node("node4", lambda state: {"my_key": "hello from node4"})
 graph.add_edge(START, "node1")
 graph.add_edge("node1", "node2")
 graph.add_edge("node2", "node3")
 graph.add_edge("node3", "node4")
 graph.add_edge("node4", END)
 return graph

def test_partial_execution_from_node2_to_node3() -> None:
 checkpointer = MemorySaver()
 graph = create_graph()
 compiled_graph = graph.compile(checkpointer=checkpointer)
 compiled_graph.update_state(
 config={
 "configurable": {
 "thread_id": "1"
 }
 },
 # The state passed into node 2 - simulating the state at
 # the end of node 1
 values={"my_key": "initial_value"},
 # Update saved state as if it came from node 1
 # Execution will resume at node 2
 as_node="node1",
 )
 result = compiled_graph.invoke(
 # Resume execution by passing None
 None,
 config={"configurable": {"thread_id": "1"}},
 # Stop after node 3 so that node 4 doesn't run
 interrupt_after="node3",
 )
 assert result["my_key"] == "hello from node3"

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/test.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Application structurePrevious](/oss/python/langgraph/application-structure)[Backward compatibilityNext](/oss/python/langgraph/backward-compatibility)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
