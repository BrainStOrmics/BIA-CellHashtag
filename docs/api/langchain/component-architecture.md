# Component architecture

[Conceptual overviews](/oss/python/concepts/products)

# Component architecture
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.LangChain’s power comes from how its components work together to create sophisticated AI applications. This page provides diagrams showcasing the relationships between different components.

## [​](#core-component-ecosystem)Core component ecosystem

The diagram below shows how LangChain’s major components connect to form complete AI applications:

### [​](#how-components-connect)How components connect

Each component layer builds on the previous ones:

{' ' * (self.list_depth - 1)}- **Input processing** – Transform raw data into structured documents

{' ' * (self.list_depth - 1)}- **Embedding & storage** – Convert text into searchable vector representations

{' ' * (self.list_depth - 1)}- **Retrieval** – Find relevant information based on user queries

{' ' * (self.list_depth - 1)}- **Generation** – Use AI models to create responses, optionally with tools

{' ' * (self.list_depth - 1)}- **Orchestration** – Coordinate everything through agents and memory systems

## [​](#component-categories)Component categories

LangChain organizes components into these main categories:

| 
 | Category | Purpose | Key Components | Use Cases
 | **[Models](/oss/python/langchain/models)** | AI reasoning and generation | Chat models, LLMs, Embedding models | Text generation, reasoning, semantic understanding
 | **[Tools](/oss/python/langchain/tools)** | External capabilities | APIs, databases, etc. | Web search, data access, computations
 | **[Agents](/oss/python/langchain/agents)** | Orchestration and reasoning | ReAct agents, tool calling agents | Nondeterministic workflows, decision making
 | **[Memory](/oss/python/langchain/short-term-memory)** | Context preservation | Message history, custom state | Conversations, stateful interactions
 | **[Retrievers](/oss/python/integrations/retrievers)** | Information access | Vector retrievers, web retrievers | RAG, knowledge base search
 | **[Document processing](/oss/python/integrations/document_loaders)** | Data ingestion | Loaders, splitters, transformers | PDF processing, web scraping
 | **[Vector Stores](/oss/python/integrations/vectorstores)** | Semantic search | Chroma, Pinecone, FAISS | Similarity search, embeddings storage

## [​](#common-patterns)Common patterns

### [​](#rag-retrieval-augmented-generation)RAG (Retrieval-Augmented generation)

### [​](#agent-with-tools)Agent with tools

### [​](#multi-agent-system)Multi-agent system

## [​](#learn-more)Learn more

{' ' * (self.list_depth - 1)}- [Creating agents](/oss/python/langchain/agents)

{' ' * (self.list_depth - 1)}- [Working with tools](/oss/python/langchain/tools)

{' ' * (self.list_depth - 1)}- [Browse integrations](/oss/python/integrations/providers/overview)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/component-architecture.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Providers and modelsPrevious](/oss/python/concepts/providers-and-models)[Memory overviewNext](/oss/python/concepts/memory)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
