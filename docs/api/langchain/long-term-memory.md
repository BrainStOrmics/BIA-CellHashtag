# Long-term memory

[Advanced usage](/oss/python/langchain/guardrails)

# Long-term memory
Copy page

Add long-term memory to LangChain agents to store and recall data across conversations and sessionsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Long-term memory lets your agent store and recall information across different conversations and sessions.
Unlike [short-term memory](/oss/python/langchain/short-term-memory), which is scoped to a single thread, long-term memory persists across threads and can be recalled at any time.
Long-term memory is built on [LangGraph stores](/oss/python/langgraph/persistence#memory-store), which save data as JSON documents organized by namespace and key.

## [​](#usage)Usage

To add long-term memory to an agent, create a store and pass it to [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent):

{' ' * (self.list_depth - 1)}- InMemoryStore
{' ' * (self.list_depth - 1)}- PostgreSQL
```
from langchain.agents import create_agent
from langchain_core.runnables import Runnable
from langgraph.store.memory import InMemoryStore

# InMemoryStore saves data to an in-memory dictionary. Use a DB-backed store in production use.
store = InMemoryStore()

agent: Runnable = create_agent(
 "claude-sonnet-4-6",
 tools=[],
 store=store,
)

```

```
pip install langgraph-checkpoint-postgres

```

```
from langchain.agents import create_agent
from langchain_core.runnables import Runnable
from langgraph.store.postgres import PostgresStore # type: ignore[import-not-found]

DB_URI = "postgresql://postgres:postgres@localhost:5442/postgres?sslmode=disable"

with PostgresStore.from_conn_string(DB_URI) as store:
 store.setup()
 agent: Runnable = create_agent(
 "claude-sonnet-4-6",
 tools=[],
 store=store,
 )

```

Tools can then read from and write to the store using the `runtime.store` parameter. See [Read long-term memory in tools](#read-long-term-memory-in-tools) and [Write long-term memory from tools](#write-long-term-memory-from-tools) for examples.
For a deeper dive into memory types (semantic, episodic, procedural) and strategies for writing memories, see the [Memory conceptual guide](/oss/python/concepts/memory#long-term-memory).

## [​](#memory-storage)Memory storage

LangGraph stores long-term memories as JSON documents in a [store](/oss/python/langgraph/persistence#memory-store).
Each memory is organized under a custom `namespace` (similar to a folder) and a distinct `key` (like a file name). Namespaces often include user or org IDs or other labels that makes it easier to organize information.
This structure enables hierarchical organization of memories. Cross-namespace searching is then supported through content filters.

{' ' * (self.list_depth - 1)}- InMemoryStore
{' ' * (self.list_depth - 1)}- PostgreSQL
```
from collections.abc import Sequence

from langgraph.store.base import IndexConfig
from langgraph.store.memory import InMemoryStore

def embed(texts: Sequence[str]) -> list[list[float]]:
 # Replace with an actual embedding function or LangChain embeddings object
 return [[1.0, 2.0] for _ in texts]

# InMemoryStore saves data to an in-memory dictionary. Use a DB-backed store in production use.
store = InMemoryStore(index=IndexConfig(embed=embed, dims=2))
user_id = "my-user"
application_context = "chitchat"
namespace = (user_id, application_context)
store.put(
 namespace,
 "a-memory",
 {
 "rules": [
 "User likes short, direct language",
 "User only speaks English & python",
 ],
 "my-key": "my-value",
 },
)
# get the "memory" by ID
item = store.get(namespace, "a-memory")
# search for "memories" within this namespace, filtering on content equivalence, sorted by vector similarity
items = store.search(
 namespace, filter={"my-key": "my-value"}, query="language preferences"
)

```

```
from collections.abc import Sequence

from langgraph.store.base import IndexConfig
from langgraph.store.postgres import PostgresStore # type: ignore[import-not-found]

def embed(texts: Sequence[str]) -> list[list[float]]:
 # Replace with an actual embedding function or LangChain embeddings object
 return [[1.0, 2.0] for _ in texts]

DB_URI = "postgresql://postgres:postgres@localhost:5442/postgres?sslmode=disable"

with PostgresStore.from_conn_string(
 DB_URI,
 index=IndexConfig(embed=embed, dims=2), # type: ignore[arg-type]
) as store:
 store.setup()
 user_id = "my-user"
 application_context = "chitchat"
 namespace = (user_id, application_context)
 store.put(
 namespace,
 "a-memory",
 {
 "rules": [
 "User likes short, direct language",
 "User only speaks English & python",
 ],
 "my-key": "my-value",
 },
 )
 item = store.get(namespace, "a-memory")
 items = store.search(
 namespace, filter={"my-key": "my-value"}, query="language preferences"
 )

```

For more information about the memory store, see the [Persistence](/oss/python/langgraph/persistence#memory-store) guide.

## [​](#read-long-term-memory-in-tools)Read long-term memory in tools

{' ' * (self.list_depth - 1)}- InMemoryStore
{' ' * (self.list_depth - 1)}- PostgreSQLGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from dataclasses import dataclass

from langchain.agents import create_agent
from langchain.tools import ToolRuntime, tool
from langchain_core.runnables import Runnable
from langgraph.store.memory import InMemoryStore

@dataclass
class Context:
 user_id: str

# InMemoryStore saves data to an in-memory dictionary. Use a DB-backed store in production.
store = InMemoryStore()

# Write sample data to the store using the put method
store.put(
 (
 "users",
 ), # Namespace to group related data together (users namespace for user data)
 "user_123", # Key within the namespace (user ID as key)
 {
 "name": "John Smith",
 "language": "English",
 }, # Data to store for the given user
)

@tool
def get_user_info(runtime: ToolRuntime[Context]) -> str:
 """Look up user info."""
 # Access the store - same as that provided to `create_agent`
 assert runtime.store is not None
 user_id = runtime.context.user_id
 # Retrieve data from store - returns StoreValue object with value and metadata
 user_info = runtime.store.get(("users",), user_id)
 return str(user_info.value) if user_info else "Unknown user"

agent: Runnable = create_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[get_user_info],
 # Pass store to agent - enables agent to access store when running tools
 store=store,
 context_schema=Context,
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "look up user information"}]},
 context=Context(user_id="user_123"),
)

```

```
from dataclasses import dataclass

from langchain.agents import create_agent
from langchain.tools import ToolRuntime, tool
from langchain_core.runnables import Runnable
from langgraph.store.postgres import PostgresStore # type: ignore[import-not-found]

@dataclass
class Context:
 user_id: str

DB_URI = "postgresql://postgres:postgres@localhost:5442/postgres?sslmode=disable"

with PostgresStore.from_conn_string(DB_URI) as store:
 store.setup()
 store.put(("users",), "user_123", {"name": "John Smith", "language": "English"})

 @tool
 def get_user_info(runtime: ToolRuntime[Context]) -> str:
 """Look up user info."""
 assert runtime.store is not None
 user_info = runtime.store.get(("users",), runtime.context.user_id)
 return str(user_info.value) if user_info else "Unknown user"

 agent: Runnable = create_agent(
 "claude-sonnet-4-6",
 tools=[get_user_info],
 store=store,
 context_schema=Context,
 )

 result = agent.invoke(
 {"messages": [{"role": "user", "content": "look up user information"}]},
 context=Context(user_id="user_123"),
 )

```

[]()

## [​](#write-long-term-memory-from-tools)Write long-term memory from tools

{' ' * (self.list_depth - 1)}- InMemoryStore
{' ' * (self.list_depth - 1)}- PostgreSQLGoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from dataclasses import dataclass

from langchain.agents import create_agent
from langchain.tools import ToolRuntime, tool
from langchain_core.runnables import Runnable
from langgraph.store.memory import InMemoryStore
from typing_extensions import TypedDict

# InMemoryStore saves data to an in-memory dictionary. Use a DB-backed store in production.
store = InMemoryStore()

@dataclass
class Context:
 user_id: str

# TypedDict defines the structure of user information for the LLM
class UserInfo(TypedDict):
 name: str

# Tool that allows agent to update user information (useful for chat applications)
@tool
def save_user_info(user_info: UserInfo, runtime: ToolRuntime[Context]) -> str:
 """Save user info."""
 # Access the store - same as that provided to `create_agent`
 assert runtime.store is not None
 store = runtime.store
 user_id = runtime.context.user_id
 # Store data in the store (namespace, key, data)
 store.put(("users",), user_id, dict(user_info))
 return "Successfully saved user info."

agent: Runnable = create_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[save_user_info],
 store=store,
 context_schema=Context,
)

# Run the agent
agent.invoke(
 {"messages": [{"role": "user", "content": "My name is John Smith"}]},
 # user_id passed in context to identify whose information is being updated
 context=Context(user_id="user_123"),
)

# You can access the store directly to get the value
item = store.get(("users",), "user_123")

```

```
from dataclasses import dataclass

from langchain.agents import create_agent
from langchain.tools import ToolRuntime, tool
from langchain_core.runnables import Runnable
from langgraph.store.postgres import PostgresStore # type: ignore[import-not-found]
from typing_extensions import TypedDict

@dataclass
class Context:
 user_id: str

class UserInfo(TypedDict):
 name: str

@tool
def save_user_info(user_info: UserInfo, runtime: ToolRuntime[Context]) -> str:
 """Save user info."""
 assert runtime.store is not None
 runtime.store.put(("users",), runtime.context.user_id, dict(user_info))
 return "Successfully saved user info."

DB_URI = "postgresql://postgres:postgres@localhost:5442/postgres?sslmode=disable"

with PostgresStore.from_conn_string(DB_URI) as store:
 store.setup()
 agent: Runnable = create_agent(
 "claude-sonnet-4-6",
 tools=[save_user_info],
 store=store,
 context_schema=Context,
 )

 agent.invoke(
 {"messages": [{"role": "user", "content": "My name is John Smith"}]},
 context=Context(user_id="user_123"),
 )

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/long-term-memory.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[RetrievalPrevious](/oss/python/langchain/retrieval)[LangSmith StudioNext](/oss/python/langchain/studio)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
