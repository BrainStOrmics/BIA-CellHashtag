# Persistence

[Capabilities](/oss/python/langgraph/persistence)

# Persistence
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.LangGraph has a built-in persistence layer that saves graph state as checkpoints. When you compile a graph with a checkpointer, a snapshot of the graph state is saved at every step of execution, organized into threads. This enables human-in-the-loop workflows, conversational memory, time travel debugging, and fault-tolerant execution.

![Checkpoints](https://mintcdn.com/langchain-5e9cc07a/-_xGPoyjhyiDWTPJ/oss/images/checkpoints.jpg?fit=max&auto=format&n=-_xGPoyjhyiDWTPJ&q=85&s=966566aaae853ed4d240c2d0d067467c)
**Agent Server handles checkpointing automatically**
When using the [Agent Server](/langsmith/agent-server), you don’t need to implement or configure checkpointers manually. The server handles all persistence infrastructure for you behind the scenes.
Trace checkpointed state and debug how your agent resumes across sessions with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-persistence). Follow the [tracing quickstart](/langsmith/trace-with-langgraph) to get set up.

## [​](#why-use-persistence)Why use persistence

Persistence is required for the following features:

{' ' * (self.list_depth - 1)}- **Human-in-the-loop**: Checkpointers facilitate [human-in-the-loop workflows](/oss/python/langgraph/interrupts) by allowing humans to inspect, interrupt, and approve graph steps. Checkpointers are needed for these workflows as the person has to be able to view the state of a graph at any point in time, and the graph has to be able to resume execution after the person has made any updates to the state. See [Interrupts](/oss/python/langgraph/interrupts) for examples.

{' ' * (self.list_depth - 1)}- **Memory**: Checkpointers allow for [“memory”](/oss/python/concepts/memory) between interactions. In the case of repeated human interactions (like conversations) any follow up messages can be sent to that thread, which will retain its memory of previous ones. See [Add memory](/oss/python/langgraph/add-memory) for information on how to add and manage conversation memory using checkpointers.

{' ' * (self.list_depth - 1)}- **Time travel**: Checkpointers allow for [“time travel”](/oss/python/langgraph/use-time-travel), allowing users to replay prior graph executions to review and / or debug specific graph steps. In addition, checkpointers make it possible to fork the graph state at arbitrary checkpoints to explore alternative trajectories.

{' ' * (self.list_depth - 1)}- **Fault-tolerance**: Checkpointing provides fault-tolerance and error recovery: if one or more nodes fail at a given superstep, you can restart your graph from the last successful step.

[]()

{' ' * (self.list_depth - 1)}- **Pending writes**: When a graph node fails mid-execution at a given [super-step](#super-steps), LangGraph stores pending checkpoint writes from any other nodes that completed successfully at that super-step. When you resume graph execution from that super-step you don’t re-run the successful nodes.

## [​](#core-concepts)Core concepts

### [​](#threads)Threads

A thread is a unique ID or thread identifier assigned to each checkpoint saved by a checkpointer. It contains the accumulated state of a sequence of [runs](/langsmith/runs). When a run is executed, the [state](/oss/python/langgraph/graph-api#state) of the underlying graph of the assistant will be persisted to the thread.
When invoking a graph with a checkpointer, you **must** specify a `thread_id` as part of the `configurable` portion of the config:

```
{"configurable": {"thread_id": "1"}}

```

A thread’s current and historical state can be retrieved. To persist state, a thread must be created prior to executing a run. The LangSmith API provides several endpoints for creating and managing threads and thread state. See the [API reference](https://reference.langchain.com/python/langsmith/) for more details.
The checkpointer uses `thread_id` as the primary key for storing and retrieving checkpoints. Without it, the checkpointer cannot save state or resume execution after an [interrupt](/oss/python/langgraph/interrupts), since the checkpointer uses `thread_id` to load the saved state.

### [​](#checkpoints)Checkpoints

The state of a thread at a particular point in time is called a checkpoint. A checkpoint is a snapshot of the graph state saved at each [super-step](#super-steps) and is represented by a `StateSnapshot` object (see [StateSnapshot fields](#statesnapshot-fields) for the full field reference).

#### [​](#super-steps)Super-steps

LangGraph creates a checkpoint at each **super-step** boundary. A super-step is a single “tick” of the graph where all nodes scheduled for that step execute (potentially in parallel). For a sequential graph like `START -> A -> B -> END`, there are separate super-steps for the input, node A, and node B — producing a checkpoint after each one. Understanding super-step boundaries is important for [time travel](/oss/python/langgraph/use-time-travel), because you can only resume execution from a checkpoint (i.e., a super-step boundary).
In addition to super-step checkpoints, LangGraph also persists writes at the **node (task) level**. As each node within a super-step finishes, its outputs are written to the checkpointer’s `checkpoint_writes` table as task entries linked to the in-progress checkpoint. These per-task writes are what enable [pending writes](#pending-writes) recovery: if another node in the same super-step fails, the successful nodes’ writes are already durable and don’t need to be re-run on resume. The full state snapshot is then committed once the super-step completes.
LangGraph also persists writes from individual node executions within a super-step. These writes are stored as tasks and used for fault tolerance: if another node in the same super-step fails, successful node writes do not need to be recomputed when you resume. These task writes are not full `StateSnapshot` checkpoints, so time travel resumes from full checkpoints at super-step boundaries.
Checkpoints are persisted and can be used to restore the state of a thread at a later time.
Let’s see what checkpoints are saved when a simple graph is invoked as follows:

```
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import InMemorySaver
from langchain_core.runnables import RunnableConfig
from typing import Annotated
from typing_extensions import TypedDict
from operator import add

class State(TypedDict):
 foo: str
 bar: Annotated[list[str], add]

def node_a(state: State):
 return {"foo": "a", "bar": ["a"]}

def node_b(state: State):
 return {"foo": "b", "bar": ["b"]}

workflow = StateGraph(State)
workflow.add_node(node_a)
workflow.add_node(node_b)
workflow.add_edge(START, "node_a")
workflow.add_edge("node_a", "node_b")
workflow.add_edge("node_b", END)

checkpointer = InMemorySaver()
graph = workflow.compile(checkpointer=checkpointer)

config: RunnableConfig = {"configurable": {"thread_id": "1"}}
graph.invoke({"foo": "", "bar":[]}, config)

```

After we run the graph, we expect to see exactly 4 checkpoints:

{' ' * (self.list_depth - 1)}- Empty checkpoint with [`START`](https://reference.langchain.com/python/langgraph/constants/START) as the next node to be executed

{' ' * (self.list_depth - 1)}- Checkpoint with the user input `{'foo': '', 'bar': []}` and `node_a` as the next node to be executed

{' ' * (self.list_depth - 1)}- Checkpoint with the outputs of `node_a` `{'foo': 'a', 'bar': ['a']}` and `node_b` as the next node to be executed

{' ' * (self.list_depth - 1)}- Checkpoint with the outputs of `node_b` `{'foo': 'b', 'bar': ['a', 'b']}` and no next nodes to be executed

Note that we `bar` channel values contain outputs from both nodes as we have a reducer for `bar` channel.

#### [​](#checkpoint-namespace)Checkpoint namespace

Each checkpoint has a `checkpoint_ns` (checkpoint namespace) field that identifies which graph or subgraph it belongs to:

{' ' * (self.list_depth - 1)}- **`""`** (empty string): The checkpoint belongs to the parent (root) graph.

{' ' * (self.list_depth - 1)}- **`"node_name:uuid"`**: The checkpoint belongs to a subgraph invoked as the given node. For nested subgraphs, namespaces are joined with `|` separators (e.g., `"outer_node:uuid|inner_node:uuid"`).

You can access the checkpoint namespace from within a node via the config:

```
from langchain_core.runnables import RunnableConfig

def my_node(state: State, config: RunnableConfig):
 checkpoint_ns = config["configurable"]["checkpoint_ns"]
 # "" for the parent graph, "node_name:uuid" for a subgraph

```

See [Subgraphs](/oss/python/langgraph/use-subgraphs) for more details on working with subgraph state and checkpoints.

## [​](#get-and-update-state)Get and update state

### [​](#get-state)Get state

When interacting with the saved graph state, you **must** specify a [thread identifier](#threads). You can view the *latest* state of the graph by calling `graph.get_state(config)`. This will return a `StateSnapshot` object that corresponds to the latest checkpoint associated with the thread ID provided in the config or a checkpoint associated with a checkpoint ID for the thread, if provided.

```
# get the latest state snapshot
config = {"configurable": {"thread_id": "1"}}
graph.get_state(config)

# get a state snapshot for a specific checkpoint_id
config = {"configurable": {"thread_id": "1", "checkpoint_id": "1ef663ba-28fe-6528-8002-5a559208592c"}}
graph.get_state(config)

```

In our example, the output of `get_state` will look like this:

```
StateSnapshot(
 values={'foo': 'b', 'bar': ['a', 'b']},
 next=(),
 config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28fe-6528-8002-5a559208592c'}},
 metadata={'source': 'loop', 'writes': {'node_b': {'foo': 'b', 'bar': ['b']}}, 'step': 2},
 created_at='2024-08-29T19:19:38.821749+00:00',
 parent_config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f9-6ec4-8001-31981c2c39f8'}}, tasks=()
)

```

#### [​](#statesnapshot-fields)StateSnapshot fields

| 
 | Field | Type | Description
 | `values` | `dict` | State channel values at this checkpoint.
 | `next` | `tuple[str, ...]` | Node names to execute next. Empty `()` means the graph is complete.
 | `config` | `dict` | Contains `thread_id`, `checkpoint_ns`, and `checkpoint_id`.
 | `metadata` | `dict` | Execution metadata. Contains `source` (`"input"`, `"loop"`, or `"update"`), `writes` (node outputs), and `step` (super-step counter).
 | `created_at` | `str` | ISO 8601 timestamp of when this checkpoint was created.
 | `parent_config` | `dict | None` | Config of the previous checkpoint. `None` for the first checkpoint.
 | `tasks` | `tuple[PregelTask, ...]` | Tasks to execute at this step. Each task has `id`, `name`, `error`, `interrupts`, and optionally `state` (subgraph snapshot, when using `subgraphs=True`).

### [​](#get-state-history)Get state history

You can get the full history of the graph execution for a given thread by calling [`graph.get_state_history(config)`](https://reference.langchain.com/python/langgraph/graphs/#langgraph.graph.state.CompiledStateGraph.get_state_history). This will return a list of `StateSnapshot` objects associated with the thread ID provided in the config. Importantly, the checkpoints will be ordered chronologically with the most recent checkpoint / `StateSnapshot` being the first in the list.

```
config = {"configurable": {"thread_id": "1"}}
list(graph.get_state_history(config))

```

In our example, the output of [`get_state_history`](https://reference.langchain.com/python/langgraph/graphs/#langgraph.graph.state.CompiledStateGraph.get_state_history) will look like this:

```
[
 StateSnapshot(
 values={'foo': 'b', 'bar': ['a', 'b']},
 next=(),
 config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28fe-6528-8002-5a559208592c'}},
 metadata={'source': 'loop', 'writes': {'node_b': {'foo': 'b', 'bar': ['b']}}, 'step': 2},
 created_at='2024-08-29T19:19:38.821749+00:00',
 parent_config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f9-6ec4-8001-31981c2c39f8'}},
 tasks=(),
 ),
 StateSnapshot(
 values={'foo': 'a', 'bar': ['a']},
 next=('node_b',),
 config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f9-6ec4-8001-31981c2c39f8'}},
 metadata={'source': 'loop', 'writes': {'node_a': {'foo': 'a', 'bar': ['a']}}, 'step': 1},
 created_at='2024-08-29T19:19:38.819946+00:00',
 parent_config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f4-6b4a-8000-ca575a13d36a'}},
 tasks=(PregelTask(id='6fb7314f-f114-5413-a1f3-d37dfe98ff44', name='node_b', error=None, interrupts=()),),
 ),
 StateSnapshot(
 values={'foo': '', 'bar': []},
 next=('node_a',),
 config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f4-6b4a-8000-ca575a13d36a'}},
 metadata={'source': 'loop', 'writes': None, 'step': 0},
 created_at='2024-08-29T19:19:38.817813+00:00',
 parent_config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f0-6c66-bfff-6723431e8481'}},
 tasks=(PregelTask(id='f1b14528-5ee5-579c-949b-23ef9bfbed58', name='node_a', error=None, interrupts=()),),
 ),
 StateSnapshot(
 values={'bar': []},
 next=('__start__',),
 config={'configurable': {'thread_id': '1', 'checkpoint_ns': '', 'checkpoint_id': '1ef663ba-28f0-6c66-bfff-6723431e8481'}},
 metadata={'source': 'input', 'writes': {'foo': ''}, 'step': -1},
 created_at='2024-08-29T19:19:38.816205+00:00',
 parent_config=None,
 tasks=(PregelTask(id='6d27aa2e-d72b-5504-a36f-8620e54a76dd', name='__start__', error=None, interrupts=()),),
 )
]

```

![State](https://mintcdn.com/langchain-5e9cc07a/-_xGPoyjhyiDWTPJ/oss/images/get_state.jpg?fit=max&auto=format&n=-_xGPoyjhyiDWTPJ&q=85&s=38ffff52be4d8806b287836295a3c058)

#### [​](#find-a-specific-checkpoint)Find a specific checkpoint

You can filter the state history to find checkpoints matching specific criteria:

```
history = list(graph.get_state_history(config))

# Find the checkpoint before a specific node executed
before_node_b = next(s for s in history if s.next == ("node_b",))

# Find a checkpoint by step number
step_2 = next(s for s in history if s.metadata["step"] == 2)

# Find checkpoints created by update_state
forks = [s for s in history if s.metadata["source"] == "update"]

# Find the checkpoint where an interrupt occurred
interrupted = next(
 s for s in history
 if s.tasks and any(t.interrupts for t in s.tasks)
)

```

### [​](#replay)Replay

Replay re-executes steps from a prior checkpoint. Invoke the graph with a prior `checkpoint_id` to re-run nodes after that checkpoint. Nodes before the checkpoint are skipped (their results are already saved). Nodes after the checkpoint re-execute, including any LLM calls, API requests, or [interrupts](/oss/python/langgraph/interrupts) — which are always re-triggered during replay.
See [Time travel](/oss/python/langgraph/use-time-travel) for full details and code examples on replaying past executions.

![Replay](https://mintcdn.com/langchain-5e9cc07a/dL5Sn6Cmy9pwtY0V/oss/images/re_play.png?fit=max&auto=format&n=dL5Sn6Cmy9pwtY0V&q=85&s=d7b34b85c106e55d181ae1f4afb50251)

### [​](#update-state)Update state

You can edit the graph state using [`update_state`](https://reference.langchain.com/python/langgraph/graphs/#langgraph.graph.state.CompiledStateGraph.update_state). This creates a new checkpoint with the updated values — it does not modify the original checkpoint. The update is treated the same as a node update: values are passed through [reducer](/oss/python/langgraph/graph-api#reducers) functions when defined, so channels with reducers *accumulate* values rather than overwrite them.
You can optionally specify `as_node` to control which node the update is treated as coming from, which affects which node executes next. See [Time travel: `as_node`](/oss/python/langgraph/use-time-travel#from-a-specific-node) for details.

![Update](https://mintcdn.com/langchain-5e9cc07a/-_xGPoyjhyiDWTPJ/oss/images/checkpoints_full_story.jpg?fit=max&auto=format&n=-_xGPoyjhyiDWTPJ&q=85&s=a52016b2c44b57bd395d6e1eac47aa36)

## [​](#memory-store)Memory store

![Model of shared state](https://mintcdn.com/langchain-5e9cc07a/dL5Sn6Cmy9pwtY0V/oss/images/shared_state.png?fit=max&auto=format&n=dL5Sn6Cmy9pwtY0V&q=85&s=354526fb48c5eb11b4b2684a2df40d6c)
A [state schema](/oss/python/langgraph/graph-api#schema) specifies a set of keys that are populated as a graph is executed. As discussed above, state can be written by a checkpointer to a thread at each graph step, enabling state persistence.
What if we want to retain some information *across threads*? Consider the case of a chatbot where we want to retain specific information about the user across *all* chat conversations (e.g., threads) with that user!
With checkpointers alone, we cannot share information across threads. This motivates the need for the [`Store`](https://reference.langchain.com/python/langgraph/store/) interface. As an illustration, we can define an `InMemoryStore` to store information about a user across threads. We simply compile our graph with a checkpointer, as before, and pass the store.
**LangGraph API handles stores automatically**
When using the LangGraph API, you don’t need to implement or configure stores manually. The API handles all storage infrastructure for you behind the scenes.
[InMemoryStore](https://reference.langchain.com/python/langchain-core/stores/InMemoryStore) is suitable for development and testing. For production, use a persistent store like `PostgresStore`, `MongoDBStore`, or `RedisStore`. All implementations extend [BaseStore](https://reference.langchain.com/python/langchain-core/stores/BaseStore), which is the type annotation to use in node function signatures.

### [​](#basic-usage)Basic usage

First, let’s showcase this in isolation without using LangGraph.

```
from langgraph.store.memory import InMemoryStore
store = InMemoryStore()

```

Memories are namespaced by a `tuple`, which in this specific example will be `(<user_id>, "memories")`. The namespace can be any length and represent anything, does not have to be user specific.

```
user_id = "1"
namespace_for_memory = (user_id, "memories")

```

We use the `store.put` method to save memories to our namespace in the store. When we do this, we specify the namespace, as defined above, and a key-value pair for the memory: the key is simply a unique identifier for the memory (`memory_id`) and the value (a dictionary) is the memory itself.

```
memory_id = str(uuid.uuid4())
memory = {"food_preference" : "I like pizza"}
store.put(namespace_for_memory, memory_id, memory)

```

We can read out memories in our namespace using the `store.search` method, which will return memories for a given user as a list, up to the `limit` argument (default `10`). With `InMemoryStore`, items are returned in insertion order, so the most recent memory is last in the list; other backends may order differently (see [Listing items in a namespace](#listing-items-in-a-namespace)).

```
memories = store.search(namespace_for_memory)
memories[-1].dict()
{'value': {'food_preference': 'I like pizza'},
 'key': '07e0caf4-1631-47b7-b15f-65515d4c1843',
 'namespace': ['1', 'memories'],
 'created_at': '2024-10-02T17:22:31.590602+00:00',
 'updated_at': '2024-10-02T17:22:31.590605+00:00'}

```

Each memory type is a Python class ([`Item`](https://langchain-ai.github.io/langgraph/reference/store/#langgraph.store.base.Item)) with certain attributes. We can access it as a dictionary by converting via `.dict` as above.
The attributes it has are:

{' ' * (self.list_depth - 1)}- 
`value`: The value (itself a dictionary) of this memory

{' ' * (self.list_depth - 1)}- 
`key`: A unique key for this memory in this namespace

{' ' * (self.list_depth - 1)}- 
`namespace`: A tuple of strings, the namespace of this memory type
While the type is `tuple[str, ...]`, it may be serialized as a list when converted to JSON (for example, `['1', 'memories']`).

{' ' * (self.list_depth - 1)}- 
`created_at`: Timestamp for when this memory was created

{' ' * (self.list_depth - 1)}- 
`updated_at`: Timestamp for when this memory was updated

### [​](#listing-items-in-a-namespace)Listing items in a namespace

Calling [`store.search`](https://reference.langchain.com/python/langgraph/store/#langgraph.store.base.BaseStore.search) (or the async [`store.asearch`](https://reference.langchain.com/python/langgraph/store/#langgraph.store.base.BaseStore.asearch)) with no `query` and no `filter` returns the items stored under `namespace_prefix`, up to `limit`. Use this to enumerate everything in a namespace when you don’t need semantic ranking.

```
# Return up to 100 items stored under ("alice", "memories").
items = store.search(("alice", "memories"), limit=100)

```

Three behaviors to keep in mind:

{' ' * (self.list_depth - 1)}- **`namespace_prefix` matches by prefix, not exactly.** `("alice",)` also returns items under `("alice", "memories")`, `("alice", "preferences")`, and so on. To restrict to a single level, pass the full namespace or filter the returned items client-side on `item.namespace`.

{' ' * (self.list_depth - 1)}- **Results past `limit` are silently truncated.** There is no overflow signal—set `limit` above your expected maximum, or paginate with `offset`.

{' ' * (self.list_depth - 1)}- **Default ordering depends on the store backend.** `PostgresStore` and `AsyncPostgresStore` return results ordered by `updated_at` descending (most recently updated first). `InMemoryStore` returns results in insertion order (most recently inserted last). Do not rely on a specific order across implementations—sort client-side on `item.updated_at` if order matters.

To page through a large namespace:

```
page_size = 50
offset = 0
while True:
 page = store.search(("alice", "memories"), limit=page_size, offset=offset)
 if not page:
 break
 for item in page:
 pass
 offset += page_size

```

To discover which namespaces exist (for example, to iterate over every user before listing their memories), use [`store.list_namespaces`](https://reference.langchain.com/python/langgraph/store/#langgraph.store.base.BaseStore.list_namespaces) or [`store.alist_namespaces`](https://reference.langchain.com/python/langgraph/store/#langgraph.store.base.BaseStore.alist_namespaces):

```
# All namespaces that start with ("alice",), truncated to two levels deep.
namespaces = store.list_namespaces(prefix=("alice",), max_depth=2)

```

### [​](#semantic-search)Semantic search

Beyond simple retrieval, the store also supports semantic search, allowing you to find memories based on meaning rather than exact matches. To enable this, configure the store with an embedding model:

```
from langchain.embeddings import init_embeddings

store = InMemoryStore(
 index={
 "embed": init_embeddings("openai:text-embedding-3-small"), # Embedding provider
 "dims": 1536, # Embedding dimensions
 "fields": ["food_preference", "$"] # Fields to embed
 }
)

```

Now when searching, you can use natural language queries to find relevant memories:

```
# Find memories about food preferences
# (This can be done after putting memories into the store)
memories = store.search(
 namespace_for_memory,
 query="What does the user like to eat?",
 limit=3 # Return top 3 matches
)

```

You can control which parts of your memories get embedded by configuring the `fields` parameter or by specifying the `index` parameter when storing memories:

```
# Store with specific fields to embed
store.put(
 namespace_for_memory,
 str(uuid.uuid4()),
 {
 "food_preference": "I love Italian cuisine",
 "context": "Discussing dinner plans"
 },
 index=["food_preference"] # Only embed "food_preferences" field
)

# Store without embedding (still retrievable, but not searchable)
store.put(
 namespace_for_memory,
 str(uuid.uuid4()),
 {"system_info": "Last updated: 2024-01-01"},
 index=False
)

```

### [​](#using-in-langgraph)Using in LangGraph

With this all in place, we use the store in LangGraph. The store works hand-in-hand with the checkpointer: the checkpointer saves state to threads, as discussed above, and the store allows us to store arbitrary information for access *across* threads. We compile the graph with both the checkpointer and the store as follows.

```
from dataclasses import dataclass
from langgraph.checkpoint.memory import InMemorySaver

@dataclass
class Context:
 user_id: str

# We need this because we want to enable threads (conversations)
checkpointer = InMemorySaver()

# ... Define the graph ...

# Compile the graph with the checkpointer and store
builder = StateGraph(MessagesState, context_schema=Context)
# ... add nodes and edges ...
graph = builder.compile(checkpointer=checkpointer, store=store)

```

We invoke the graph with a `thread_id`, as before, and also with a `user_id`, which we’ll use to namespace our memories to this particular user as we showed above.

```
# Invoke the graph
config = {"configurable": {"thread_id": "1"}}

# First let's just say hi to the AI
for update in graph.stream(
 {"messages": [{"role": "user", "content": "hi"}]},
 config,
 stream_mode="updates",
 context=Context(user_id="1"),
):
 print(update)

```

You can access the store and the `user_id` in *any node* by using the `Runtime` object. The `Runtime` is automatically injected by LangGraph when you add it as a parameter to your node function. Here’s how you might use it to save memories:

```
from langgraph.runtime import Runtime
from dataclasses import dataclass

@dataclass
class Context:
 user_id: str

async def update_memory(state: MessagesState, runtime: Runtime[Context]):

 # Get the user id from the runtime context
 user_id = runtime.context.user_id

 # Namespace the memory
 namespace = (user_id, "memories")

 # ... Analyze conversation and create a new memory

 # Create a new memory ID
 memory_id = str(uuid.uuid4())

 # We create a new memory
 await runtime.store.aput(namespace, memory_id, {"memory": memory})

```

As we showed above, we can also access the store in any node and use the `store.search` method to get memories. Recall the memories are returned as a list of objects that can be converted to a dictionary.

```
memories[-1].dict()
{'value': {'food_preference': 'I like pizza'},
 'key': '07e0caf4-1631-47b7-b15f-65515d4c1843',
 'namespace': ['1', 'memories'],
 'created_at': '2024-10-02T17:22:31.590602+00:00',
 'updated_at': '2024-10-02T17:22:31.590605+00:00'}

```

We can access the memories and use them in our model call.

```
from dataclasses import dataclass
from langgraph.runtime import Runtime

@dataclass
class Context:
 user_id: str

async def call_model(state: MessagesState, runtime: Runtime[Context]):
 # Get the user id from the runtime context
 user_id = runtime.context.user_id

 # Namespace the memory
 namespace = (user_id, "memories")

 # Search based on the most recent message
 memories = await runtime.store.asearch(
 namespace,
 query=state["messages"][-1].content,
 limit=3
 )
 info = "\n".join([d.value["memory"] for d in memories])

 # ... Use memories in the model call

```

If we create a new thread, we can still access the same memories so long as the `user_id` is the same.

```
# Invoke the graph on a new thread
config = {"configurable": {"thread_id": "2"}}

# Let's say hi again
for update in graph.stream(
 {"messages": [{"role": "user", "content": "hi, tell me about my memories"}]},
 config,
 stream_mode="updates",
 context=Context(user_id="1"),
):
 print(update)

```

When we use the LangSmith, either locally (e.g., in [Studio](/langsmith/studio)) or [hosted with LangSmith](/langsmith/platform-setup), the base store is available to use by default and does not need to be specified during graph compilation. To enable semantic search, however, you **do** need to configure the indexing settings in your `langgraph.json` file. For example:

```
{
 ...
 "store": {
 "index": {
 "embed": "openai:text-embeddings-3-small",
 "dims": 1536,
 "fields": ["$"]
 }
 }
}

```

See the [deployment guide](/langsmith/semantic-search) for more details and configuration options.

## [​](#optimize-checkpoint-storage)Optimize checkpoint storage

By default, LangGraph checkpoints write the full value of every state channel at each super-step. For long-running threads with large accumulations—such as multi-turn conversations—this can produce significant storage growth over time.
[`DeltaChannel`](https://reference.langchain.com/python/langgraph/channels/delta/DeltaChannel) stores only incremental deltas instead of the full accumulated value, substantially reducing checkpoint size for append-heavy channels. See [DeltaChannel](/oss/python/langgraph/pregel#deltachannel-beta) for usage and the storage-vs-latency tradeoff.
`DeltaChannel` requires `langgraph>=1.2` and is currently in beta. The API may change in future releases.

## [​](#checkpointer-libraries)Checkpointer libraries

Under the hood, checkpointing is powered by checkpointer objects that conform to [`BaseCheckpointSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.base.BaseCheckpointSaver) interface. LangGraph provides several checkpointer implementations, all implemented via standalone, installable libraries.
See [checkpointer integrations](/oss/python/integrations/checkpointers/index) for available providers.

{' ' * (self.list_depth - 1)}- `langgraph-checkpoint`: The base interface for checkpointer savers ([`BaseCheckpointSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.base.BaseCheckpointSaver)) and serialization/deserialization interface ([`SerializerProtocol`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.base.SerializerProtocol)). Includes in-memory checkpointer implementation ([`InMemorySaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.memory.InMemorySaver)) for experimentation. LangGraph comes with `langgraph-checkpoint` included.

{' ' * (self.list_depth - 1)}- `langgraph-checkpoint-sqlite`: An implementation of LangGraph checkpointer that uses SQLite database ([`SqliteSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.sqlite.SqliteSaver) / [`AsyncSqliteSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.sqlite.aio.AsyncSqliteSaver)). Ideal for experimentation and local workflows. Needs to be installed separately.

{' ' * (self.list_depth - 1)}- `langgraph-checkpoint-postgres`: An advanced checkpointer that uses Postgres database ([`PostgresSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.postgres.PostgresSaver) / [`AsyncPostgresSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.postgres.aio.AsyncPostgresSaver)), used in LangSmith. Ideal for using in production. Needs to be installed separately.

{' ' * (self.list_depth - 1)}- `langchain-azure-cosmosdb`: An implementation of LangGraph checkpointer that uses Azure Cosmos DB for NoSQL ([`CosmosDBSaverSync`](https://reference.langchain.com/python/langchain-azure-cosmosdb/) / [`CosmosDBSaver`](https://reference.langchain.com/python/langchain-azure-cosmosdb/)). Ideal for using in production with Azure. Supports both sync and async operations, with Microsoft Entra ID authentication. Needs to be installed separately.

### [​](#checkpointer-interface)Checkpointer interface

Each checkpointer conforms to [`BaseCheckpointSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.base.BaseCheckpointSaver) interface and implements the following methods:

{' ' * (self.list_depth - 1)}- `.put` - Store a checkpoint with its configuration and metadata.

{' ' * (self.list_depth - 1)}- `.put_writes` - Store intermediate writes linked to a checkpoint (i.e. [pending writes](#pending-writes)).

{' ' * (self.list_depth - 1)}- `.get_tuple` - Fetch a checkpoint tuple using for a given configuration (`thread_id` and `checkpoint_id`). This is used to populate `StateSnapshot` in `graph.get_state()`.

{' ' * (self.list_depth - 1)}- `.list` - List checkpoints that match a given configuration and filter criteria. This is used to populate state history in `graph.get_state_history()`

If the checkpointer is used with asynchronous graph execution (i.e. executing the graph via `.ainvoke`, `.astream`, `.abatch`), asynchronous versions of the above methods will be used (`.aput`, `.aput_writes`, `.aget_tuple`, `.alist`).
For running your graph asynchronously, you can use [`InMemorySaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.memory.InMemorySaver), or async versions of Sqlite/Postgres checkpointers — [`AsyncSqliteSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.sqlite.aio.AsyncSqliteSaver) / [`AsyncPostgresSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.postgres.aio.AsyncPostgresSaver) checkpointers.

### [​](#serializer)Serializer

When checkpointers save the graph state, they need to serialize the channel values in the state. This is done using serializer objects.
`langgraph_checkpoint` defines [protocol](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.base.SerializerProtocol) for implementing serializers provides a default implementation ([`JsonPlusSerializer`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.jsonplus.JsonPlusSerializer)) that handles a wide variety of types, including LangChain and LangGraph primitives, datetimes, enums and more.

#### [​](#serialization-with-pickle)Serialization with `pickle`

The default serializer, [`JsonPlusSerializer`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.jsonplus.JsonPlusSerializer), uses ormsgpack and JSON under the hood, which is not suitable for all types of objects.
If you want to fallback to pickle for objects not currently supported by our msgpack encoder (such as Pandas dataframes),
you can use the `pickle_fallback` argument of the [`JsonPlusSerializer`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.jsonplus.JsonPlusSerializer):

```
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.checkpoint.serde.jsonplus import JsonPlusSerializer

# ... Define the graph ...
graph.compile(
 checkpointer=InMemorySaver(serde=JsonPlusSerializer(pickle_fallback=True))
)

```

#### [​](#encryption)Encryption

Checkpointers can optionally encrypt all persisted state. To enable this, pass an instance of [`EncryptedSerializer`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.encrypted.EncryptedSerializer) to the `serde` argument of any [`BaseCheckpointSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.base.BaseCheckpointSaver) implementation. The easiest way to create an encrypted serializer is via [`from_pycryptodome_aes`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.encrypted.EncryptedSerializer.from_pycryptodome_aes), which reads the AES key from the `LANGGRAPH_AES_KEY` environment variable (or accepts a `key` argument):

```
import sqlite3

from langgraph.checkpoint.serde.encrypted import EncryptedSerializer
from langgraph.checkpoint.sqlite import SqliteSaver

serde = EncryptedSerializer.from_pycryptodome_aes() # reads LANGGRAPH_AES_KEY
checkpointer = SqliteSaver(sqlite3.connect("checkpoint.db"), serde=serde)

```

```
from langgraph.checkpoint.serde.encrypted import EncryptedSerializer
from langgraph.checkpoint.postgres import PostgresSaver

serde = EncryptedSerializer.from_pycryptodome_aes()
checkpointer = PostgresSaver.from_conn_string("postgresql://...", serde=serde)
checkpointer.setup()

```

When running on LangSmith, encryption is automatically enabled whenever `LANGGRAPH_AES_KEY` is present, so you only need to provide the environment variable. Other encryption schemes can be used by implementing [`CipherProtocol`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.base.CipherProtocol) and supplying it to [`EncryptedSerializer`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.serde.encrypted.EncryptedSerializer).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/persistence.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Workflows and agentsPrevious](/oss/python/langgraph/workflows-agents)[Durable executionNext](/oss/python/langgraph/durable-execution)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
