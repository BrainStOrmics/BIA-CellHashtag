# Build a custom RAG agent with LangGraph

[Tutorials](/oss/python/deepagents/data-analysis)[LangGraph](/oss/python/langgraph/agentic-rag)

# Build a custom RAG agent with LangGraph
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

In this tutorial we will build a [retrieval](/oss/python/langchain/retrieval) agent using LangGraph.
LangChain offers built-in [agent](/oss/python/langchain/agents) implementations, implemented using [LangGraph](/oss/python/langgraph/overview) primitives. If deeper customization is required, agents can be implemented directly in LangGraph. This guide demonstrates an example implementation of a retrieval agent. [Retrieval](/oss/python/langchain/retrieval) agents are useful when you want an LLM to make a decision about whether to retrieve context from a vectorstore or respond to the user directly.
By the end of the tutorial we will have done the following:

{' ' * (self.list_depth - 1)}- Fetch and preprocess documents that will be used for retrieval.

{' ' * (self.list_depth - 1)}- Index those documents for semantic search and create a retriever tool for the agent.

{' ' * (self.list_depth - 1)}- Build an agentic RAG system that can decide when to use the retriever tool.

![Hybrid RAG](https://mintcdn.com/langchain-5e9cc07a/I6RpA28iE233vhYX/images/langgraph-hybrid-rag-tutorial.png?fit=max&auto=format&n=I6RpA28iE233vhYX&q=85&s=855348219691485642b22a1419939ea7)

### [​](#concepts)Concepts

We will cover the following concepts:

{' ' * (self.list_depth - 1)}- [Retrieval](/oss/python/langchain/retrieval) using [document loaders](/oss/python/integrations/document_loaders), [text splitters](/oss/python/integrations/splitters), [embeddings](/oss/python/integrations/embeddings), and [vector stores](/oss/python/integrations/vectorstores)

{' ' * (self.list_depth - 1)}- The LangGraph [Graph API](/oss/python/langgraph/graph-api), including state, nodes, edges, and conditional edges.

## [​](#setup)Setup

Let’s download the required packages and set our API keys:

```
pip install -U langgraph "langchain[openai]" langchain-text-splitters bs4 requests

```

```
import getpass
import os

def _set_env(key: str):
 if key not in os.environ:
 os.environ[key] = getpass.getpass(f"{key}:")

_set_env("OPENAI_API_KEY")

```

Sign up for LangSmith to quickly spot issues and improve the performance of your LangGraph projects. [LangSmith](https://docs.smith.langchain.com) lets you use trace data to debug, test, and monitor your LLM apps built with LangGraph.

## [​](#1-preprocess-documents)1. Preprocess documents

{' ' * (self.list_depth - 1)}- Fetch documents to use in our RAG system. We will use three of the most recent pages from [Lilian Weng’s excellent blog](https://lilianweng.github.io/). We’ll start by fetching the content of the pages with a minimal helper built on `requests` and `BeautifulSoup`:

```
import bs4
import requests
from langchain_core.documents import Document

# Below is a minimal helper for demonstration purposes.
def load_web_page(url: str, bs_kwargs: dict | None = None) -> list[Document]:
 response = requests.get(url)
 response.raise_for_status()
 soup = bs4.BeautifulSoup(response.text, "html.parser", **(bs_kwargs or {}))
 return [Document(page_content=soup.get_text(), metadata={"source": url})]

urls = [
 "https://lilianweng.github.io/posts/2024-11-28-reward-hacking/",
 "https://lilianweng.github.io/posts/2024-07-07-hallucination/",
 "https://lilianweng.github.io/posts/2024-04-12-diffusion-video/",
]

docs = [load_web_page(url) for url in urls]

```

```
docs[0][0].page_content.strip()[:1000]

```

{' ' * (self.list_depth - 1)}- Split the fetched documents into smaller chunks for indexing into our vectorstore:

```
from langchain_text_splitters import RecursiveCharacterTextSplitter

docs_list = [item for sublist in docs for item in sublist]

text_splitter = RecursiveCharacterTextSplitter.from_tiktoken_encoder(
 chunk_size=100, chunk_overlap=50
)
doc_splits = text_splitter.split_documents(docs_list)

```

```
doc_splits[0].page_content.strip()

```

## [​](#2-create-a-retriever-tool)2. Create a retriever tool

Now that we have our split documents, we can index them into a vector store that we’ll use for semantic search.

{' ' * (self.list_depth - 1)}- Use an in-memory vector store and OpenAI embeddings:

```
from langchain_core.vectorstores import InMemoryVectorStore
from langchain_openai import OpenAIEmbeddings

vectorstore = InMemoryVectorStore.from_documents(
 documents=doc_splits, embedding=OpenAIEmbeddings()
)
retriever = vectorstore.as_retriever()

```

{' ' * (self.list_depth - 1)}- Create a retriever tool using the `@tool` decorator:

```
from langchain.tools import tool

@tool
def retrieve_blog_posts(query: str) -> str:
 """Search and return information about Lilian Weng blog posts."""
 docs = retriever.invoke(query)
 return "\n\n".join([doc.page_content for doc in docs])

retriever_tool = retrieve_blog_posts

```

{' ' * (self.list_depth - 1)}- Test the tool:

```
retriever_tool.invoke({"query": "types of reward hacking"})

```

## [​](#3-generate-query)3. Generate query

Now we will start building components ([nodes](/oss/python/langgraph/graph-api#nodes) and [edges](/oss/python/langgraph/graph-api#edges)) for our agentic RAG graph.
Note that the components will operate on the [`MessagesState`](/oss/python/langgraph/graph-api#messagesstate)—graph state that contains a `messages` key with a list of [chat messages](https://python.langchain.com/docs/concepts/messages/).

{' ' * (self.list_depth - 1)}- Build a `generate_query_or_respond` node. It will call an LLM to generate a response based on the current graph state (list of messages). Given the input messages, it will decide to retrieve using the retriever tool, or respond directly to the user. Note that we’re giving the chat model access to the `retriever_tool` we created earlier via `.bind_tools`:

```
from langgraph.graph import MessagesState
from langchain.chat_models import init_chat_model

response_model = init_chat_model("gpt-5.4", temperature=0)

def generate_query_or_respond(state: MessagesState):
 """Call the model to generate a response based on the current state. Given
 the question, it will decide to retrieve using the retriever tool, or simply respond to the user.
 """
 response = (
 response_model
 .bind_tools([retriever_tool]).invoke(state["messages"])
 )
 return {"messages": [response]}

```

{' ' * (self.list_depth - 1)}- Try it on a random input:

```
input = {"messages": [{"role": "user", "content": "hello!"}]}
generate_query_or_respond(input)["messages"][-1].pretty_print()

```

**Output:**

```
================================== Ai Message ==================================

Hello! How can I help you today?

```

{' ' * (self.list_depth - 1)}- Ask a question that requires semantic search:

```
input = {
 "messages": [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 }
 ]
}
generate_query_or_respond(input)["messages"][-1].pretty_print()

```

**Output:**

```
================================== Ai Message ==================================
Tool Calls:
retrieve_blog_posts (call_tYQxgfIlnQUDMdtAhdbXNwIM)
Call ID: call_tYQxgfIlnQUDMdtAhdbXNwIM
Args:
 query: types of reward hacking

```

## [​](#4-grade-documents)4. Grade documents

{' ' * (self.list_depth - 1)}- Add a [conditional edge](/oss/python/langgraph/graph-api#conditional-edges)—`grade_documents`—to determine whether the retrieved documents are relevant to the question. We will use a model with a structured output schema `GradeDocuments` for document grading. The `grade_documents` function will return the name of the node to go to based on the grading decision (`generate_answer` or `rewrite_question`):

```
from pydantic import BaseModel, Field
from typing import Literal

GRADE_PROMPT = (
 "You are a grader assessing relevance of a retrieved document to a user question. \n "
 "Here is the retrieved document: \n\n {context} \n\n"
 "Here is the user question: {question} \n"
 "If the document contains keyword(s) or semantic meaning related to the user question, grade it as relevant. \n"
 "Give a binary score 'yes' or 'no' score to indicate whether the document is relevant to the question."
)

class GradeDocuments(BaseModel):
 """Grade documents using a binary score for relevance check."""

 binary_score: str = Field(
 description="Relevance score: 'yes' if relevant, or 'no' if not relevant"
 )

grader_model = init_chat_model("gpt-5.4", temperature=0)

def grade_documents(
 state: MessagesState,
) -> Literal["generate_answer", "rewrite_question"]:
 """Determine whether the retrieved documents are relevant to the question."""
 question = state["messages"][0].content
 context = state["messages"][-1].content

 prompt = GRADE_PROMPT.format(question=question, context=context)
 response = (
 grader_model
 .with_structured_output(GradeDocuments).invoke(
 [{"role": "user", "content": prompt}]
 )
 )
 score = response.binary_score

 if score == "yes":
 return "generate_answer"
 else:
 return "rewrite_question"

```

{' ' * (self.list_depth - 1)}- Run this with irrelevant documents in the tool response:

```
from langchain_core.messages import convert_to_messages

input = {
 "messages": convert_to_messages(
 [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 },
 {
 "role": "assistant",
 "content": "",
 "tool_calls": [
 {
 "id": "1",
 "name": "retrieve_blog_posts",
 "args": {"query": "types of reward hacking"},
 }
 ],
 },
 {"role": "tool", "content": "meow", "tool_call_id": "1"},
 ]
 )
}
grade_documents(input)

```

{' ' * (self.list_depth - 1)}- Confirm that the relevant documents are classified as such:

```
input = {
 "messages": convert_to_messages(
 [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 },
 {
 "role": "assistant",
 "content": "",
 "tool_calls": [
 {
 "id": "1",
 "name": "retrieve_blog_posts",
 "args": {"query": "types of reward hacking"},
 }
 ],
 },
 {
 "role": "tool",
 "content": "reward hacking can be categorized into two types: environment or goal misspecification, and reward tampering",
 "tool_call_id": "1",
 },
 ]
 )
}
grade_documents(input)

```

## [​](#5-rewrite-question)5. Rewrite question

{' ' * (self.list_depth - 1)}- Build the `rewrite_question` node. The retriever tool can return potentially irrelevant documents, which indicates a need to improve the original user question. To do so, we will call the `rewrite_question` node:

```
from langchain.messages import HumanMessage

REWRITE_PROMPT = (
 "Look at the input and try to reason about the underlying semantic intent / meaning.\n"
 "Here is the initial question:"
 "\n ------- \n"
 "{question}"
 "\n ------- \n"
 "Formulate an improved question:"
)

def rewrite_question(state: MessagesState):
 """Rewrite the original user question."""
 messages = state["messages"]
 question = messages[0].content
 prompt = REWRITE_PROMPT.format(question=question)
 response = response_model.invoke([{"role": "user", "content": prompt}])
 return {"messages": [HumanMessage(content=response.content)]}

```

{' ' * (self.list_depth - 1)}- Try it out:

```
input = {
 "messages": convert_to_messages(
 [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 },
 {
 "role": "assistant",
 "content": "",
 "tool_calls": [
 {
 "id": "1",
 "name": "retrieve_blog_posts",
 "args": {"query": "types of reward hacking"},
 }
 ],
 },
 {"role": "tool", "content": "meow", "tool_call_id": "1"},
 ]
 )
}

response = rewrite_question(input)
print(response["messages"][-1].content)

```

**Output:**

```
What are the different types of reward hacking described by Lilian Weng, and how does she explain them?

```

## [​](#6-generate-an-answer)6. Generate an answer

{' ' * (self.list_depth - 1)}- Build `generate_answer` node: if we pass the grader checks, we can generate the final answer based on the original question and the retrieved context:

```
GENERATE_PROMPT = (
 "You are an assistant for question-answering tasks. "
 "Use the following pieces of retrieved context to answer the question. "
 "If you don't know the answer, just say that you don't know. "
 "Use three sentences maximum and keep the answer concise.\n"
 "Question: {question} \n"
 "Context: {context}"
)

def generate_answer(state: MessagesState):
 """Generate an answer."""
 question = state["messages"][0].content
 context = state["messages"][-1].content
 prompt = GENERATE_PROMPT.format(question=question, context=context)
 response = response_model.invoke([{"role": "user", "content": prompt}])
 return {"messages": [response]}

```

{' ' * (self.list_depth - 1)}- Try it:

```
input = {
 "messages": convert_to_messages(
 [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 },
 {
 "role": "assistant",
 "content": "",
 "tool_calls": [
 {
 "id": "1",
 "name": "retrieve_blog_posts",
 "args": {"query": "types of reward hacking"},
 }
 ],
 },
 {
 "role": "tool",
 "content": "reward hacking can be categorized into two types: environment or goal misspecification, and reward tampering",
 "tool_call_id": "1",
 },
 ]
 )
}

response = generate_answer(input)
response["messages"][-1].pretty_print()

```

**Output:**

```
================================== Ai Message ==================================

Lilian Weng categorizes reward hacking into two types: environment or goal misspecification, and reward tampering. She considers reward hacking as a broad concept that includes both of these categories. Reward hacking occurs when an agent exploits flaws or ambiguities in the reward function to achieve high rewards without performing the intended behaviors.

```

## [​](#7-assemble-the-graph)7. Assemble the graph

Now we’ll assemble all the nodes and edges into a complete graph:

{' ' * (self.list_depth - 1)}- Start with a `generate_query_or_respond` and determine if we need to call `retriever_tool`

{' ' * (self.list_depth - 1)}- Route to next step based on whether the model made tool calls:

{' ' * (self.list_depth - 1)}- If `generate_query_or_respond` returned `tool_calls`, call `retriever_tool` to retrieve context

{' ' * (self.list_depth - 1)}- Otherwise, respond directly to the user

{' ' * (self.list_depth - 1)}- Grade retrieved document content for relevance to the question (`grade_documents`) and route to next step:

{' ' * (self.list_depth - 1)}- If not relevant, rewrite the question using `rewrite_question` and then call `generate_query_or_respond` again

{' ' * (self.list_depth - 1)}- If relevant, proceed to `generate_answer` and generate final response using the [`ToolMessage`](https://reference.langchain.com/python/langchain-core/messages/tool/ToolMessage) with the retrieved document context

```
from langgraph.graph import END, START, StateGraph
from langgraph.prebuilt import ToolNode

workflow = StateGraph(MessagesState)

# Define the nodes we will cycle between
workflow.add_node(generate_query_or_respond)
workflow.add_node("retrieve", ToolNode([retriever_tool]))
workflow.add_node(rewrite_question)
workflow.add_node(generate_answer)

workflow.add_edge(START, "generate_query_or_respond")

# Route based on whether the model requested tool calls.
def route_on_tool_calls(state: MessagesState):
 last_message = state["messages"][-1]
 if getattr(last_message, "tool_calls", None):
 return "tools"
 return END

# Decide whether to retrieve
workflow.add_conditional_edges(
 "generate_query_or_respond",
 # Assess LLM decision (call `retriever_tool` tool or respond to the user)
 route_on_tool_calls,
 {
 # Translate the condition outputs to nodes in our graph
 "tools": "retrieve",
 END: END,
 },
)

# Edges taken after the `action` node is called.
workflow.add_conditional_edges(
 "retrieve",
 # Assess agent decision
 grade_documents,
)
workflow.add_edge("generate_answer", END)
workflow.add_edge("rewrite_question", "generate_query_or_respond")

# Compile
graph = workflow.compile()

```

Visualize the graph:

```
from IPython.display import Image, display

display(Image(graph.get_graph().draw_mermaid_png()))

```

![SQL agent graph](https://mintcdn.com/langchain-5e9cc07a/-_xGPoyjhyiDWTPJ/oss/images/agentic-rag-output.png?fit=max&auto=format&n=-_xGPoyjhyiDWTPJ&q=85&s=ddedbd57514888e614ece260092201df)

## [​](#8-run-the-agentic-rag)8. Run the agentic RAG

Now let’s test the complete graph by running it with a question:

```
for chunk in graph.stream(
 {
 "messages": [
 {
 "role": "user",
 "content": "What does Lilian Weng say about types of reward hacking?",
 }
 ]
 }
):
 for node, update in chunk.items():
 print("Update from node", node)
 update["messages"][-1].pretty_print()
 print("\n\n")

```

**Output:**

```
Update from node generate_query_or_respond
================================== Ai Message ==================================
Tool Calls:
 retrieve_blog_posts (call_NYu2vq4km9nNNEFqJwefWKu1)
 Call ID: call_NYu2vq4km9nNNEFqJwefWKu1
 Args:
 query: types of reward hacking

Update from node retrieve
================================= Tool Message ==================================
Name: retrieve_blog_posts

(Note: Some work defines reward tampering as a distinct category of misalignment behavior from reward hacking. But I consider reward hacking as a broader concept here.)
At a high level, reward hacking can be categorized into two types: environment or goal misspecification, and reward tampering.

Why does Reward Hacking Exist?#

Pan et al. (2022) investigated reward hacking as a function of agent capabilities, including (1) model size, (2) action space resolution, (3) observation space noise, and (4) training time. They also proposed a taxonomy of three types of misspecified proxy rewards:

Let's Define Reward Hacking#
Reward shaping in RL is challenging. Reward hacking occurs when an RL agent exploits flaws or ambiguities in the reward function to obtain high rewards without genuinely learning the intended behaviors or completing the task as designed. In recent years, several related concepts have been proposed, all referring to some form of reward hacking:

Update from node generate_answer
================================== Ai Message ==================================

Lilian Weng categorizes reward hacking into two types: environment or goal misspecification, and reward tampering. She considers reward hacking as a broad concept that includes both of these categories. Reward hacking occurs when an agent exploits flaws or ambiguities in the reward function to achieve high rewards without performing the intended behaviors.

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/agentic-rag.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a SQL assistant with on-demand skillsPrevious](/oss/python/langchain/multi-agent/skills-sql-assistant)[Build a custom SQL agentNext](/oss/python/langgraph/sql-agent)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
