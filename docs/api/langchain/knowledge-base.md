# Build a semantic search engine with LangChain

[Tutorials](/oss/python/deepagents/data-analysis)[LangChain](/oss/python/langchain/deep-agent-from-scratch)

# Build a semantic search engine with LangChain
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

This tutorial will familiarize you with LangChain’s [embedding](/oss/python/integrations/embeddings) and [vector store](/oss/python/integrations/vectorstores) abstractions. These abstractions are designed to support retrieval of data— from (vector) databases and other sources — for integration with LLM workflows. They are important for applications that fetch data to be reasoned over as part of model inference, as in the case of retrieval-augmented generation, or [RAG](/oss/python/langchain/retrieval).
Here we will build a search engine over a PDF document. This will allow us to retrieve passages in the PDF that are similar to an input query. The guide also includes a minimal RAG implementation on top of the search engine.

### [​](#concepts)Concepts

This guide focuses on retrieval of text data. We will cover the following concepts:

{' ' * (self.list_depth - 1)}- [Documents](https://reference.langchain.com/python/langchain-core/documents);

{' ' * (self.list_depth - 1)}- [Text splitters](/oss/python/integrations/splitters);

{' ' * (self.list_depth - 1)}- [Embeddings](/oss/python/integrations/embeddings);

{' ' * (self.list_depth - 1)}- [Vector stores](/oss/python/integrations/vectorstores) and [retrievers](/oss/python/integrations/retrievers).

## [​](#setup)Setup

### [​](#installation)Installation

This tutorial reads a PDF using the `pypdf` package:
pipcondauv
```
pip install pypdf

```

For more details, see our [Installation guide](/oss/python/langchain/install).

### [​](#langsmith)LangSmith

Many of the applications you build with LangChain will contain multiple steps with multiple invocations of LLM calls.
As these applications get more and more complex, it becomes crucial to be able to inspect what exactly is going on inside your chain or agent.
The best way to do this is with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-knowledge-base).
After you sign up at the link above, make sure to set your environment variables to start logging traces:

```
export LANGSMITH_TRACING="true"
export LANGSMITH_API_KEY="..."

```

Or, if in a notebook, you can set them with:

```
import getpass
import os

os.environ["LANGSMITH_TRACING"] = "true"
os.environ["LANGSMITH_API_KEY"] = getpass.getpass()

```

## [​](#1-documents)1. Documents

LangChain implements a [Document](https://reference.langchain.com/python/langchain-core/documents/base/Document) abstraction, which is intended to represent a unit of text and associated metadata. It has three attributes:

{' ' * (self.list_depth - 1)}- `page_content`: a string representing the content;

{' ' * (self.list_depth - 1)}- `metadata`: a dict containing arbitrary metadata;

{' ' * (self.list_depth - 1)}- `id`: (optional) a string identifier for the document.

The `metadata` attribute can capture information about the source of the document, its relationship to other documents, and other information. Note that an individual [`Document`](https://reference.langchain.com/python/langchain-core/documents/base/Document) object often represents a chunk of a larger document.
We can generate sample documents when desired:

```
from langchain_core.documents import Document

documents = [
 Document(
 page_content="Dogs are great companions, known for their loyalty and friendliness.",
 metadata={"source": "mammal-pets-doc"},
 ),
 Document(
 page_content="Cats are independent pets that often enjoy their own space.",
 metadata={"source": "mammal-pets-doc"},
 ),
]

```

## [​](#2-embeddings)2. Embeddings

Vector search is a common way to store and search over unstructured data (such as unstructured text). The idea is to store numeric vectors that are associated with the text. Given a query, we can [embed](/oss/python/integrations/embeddings) it as a vector of the same dimension and use vector similarity metrics (such as cosine similarity) to identify related text.
LangChain supports embeddings from [dozens of providers](/oss/python/integrations/embeddings). These models specify how text should be converted into a numeric vector. Let’s select a model:

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Azure
{' ' * (self.list_depth - 1)}- Google Gemini
{' ' * (self.list_depth - 1)}- Google Vertex
{' ' * (self.list_depth - 1)}- AWS
{' ' * (self.list_depth - 1)}- HuggingFace
{' ' * (self.list_depth - 1)}- Ollama
{' ' * (self.list_depth - 1)}- Cohere
{' ' * (self.list_depth - 1)}- MistralAI
{' ' * (self.list_depth - 1)}- Nomic
{' ' * (self.list_depth - 1)}- NVIDIA
{' ' * (self.list_depth - 1)}- Voyage AI
{' ' * (self.list_depth - 1)}- IBM watsonx
{' ' * (self.list_depth - 1)}- Fake
{' ' * (self.list_depth - 1)}- Isaacus
```
pip install -U "langchain-openai"

```

```
import getpass
import os

if not os.environ.get("OPENAI_API_KEY"):
 os.environ["OPENAI_API_KEY"] = getpass.getpass("Enter API key for OpenAI: ")

from langchain_openai import OpenAIEmbeddings

embeddings = OpenAIEmbeddings(model="text-embedding-3-large")

```

```
pip install -U "langchain-openai"

```

```
import getpass
import os

if not os.environ.get("AZURE_OPENAI_API_KEY"):
 os.environ["AZURE_OPENAI_API_KEY"] = getpass.getpass("Enter API key for Azure: ")

from langchain_openai import AzureOpenAIEmbeddings

embeddings = AzureOpenAIEmbeddings(
 azure_endpoint=os.environ["AZURE_OPENAI_ENDPOINT"],
 azure_deployment=os.environ["AZURE_OPENAI_DEPLOYMENT_NAME"],
 openai_api_version=os.environ["AZURE_OPENAI_API_VERSION"],
)

```

```
pip install -qU langchain-google-genai

```

```
import getpass
import os

if not os.environ.get("GOOGLE_API_KEY"):
 os.environ["GOOGLE_API_KEY"] = getpass.getpass("Enter API key for Google Gemini: ")

from langchain_google_genai import GoogleGenerativeAIEmbeddings

embeddings = GoogleGenerativeAIEmbeddings(model="models/gemini-embedding-001")

```

```
pip install -qU langchain-google-vertexai

```

```
from langchain_google_vertexai import VertexAIEmbeddings

embeddings = VertexAIEmbeddings(model="text-embedding-005")

```

```
pip install -qU langchain-aws

```

```
from langchain_aws import BedrockEmbeddings

embeddings = BedrockEmbeddings(model_id="amazon.titan-embed-text-v2:0")

```

```
pip install -qU langchain-huggingface

```

```
from langchain_huggingface import HuggingFaceEmbeddings

embeddings = HuggingFaceEmbeddings(
 model_name="sentence-transformers/all-mpnet-base-v2",
 encode_kwargs={"normalize_embeddings": True},
)

```

```
pip install -qU langchain-ollama

```

```
from langchain_ollama import OllamaEmbeddings

embeddings = OllamaEmbeddings(model="llama3")

```

```
pip install -qU langchain-cohere

```

```
import getpass
import os

if not os.environ.get("COHERE_API_KEY"):
 os.environ["COHERE_API_KEY"] = getpass.getpass("Enter API key for Cohere: ")

from langchain_cohere import CohereEmbeddings

embeddings = CohereEmbeddings(model="embed-english-v3.0")

```

```
pip install -qU langchain-mistralai

```

```
import getpass
import os

if not os.environ.get("MISTRALAI_API_KEY"):
 os.environ["MISTRALAI_API_KEY"] = getpass.getpass("Enter API key for MistralAI: ")

from langchain_mistralai import MistralAIEmbeddings

embeddings = MistralAIEmbeddings(model="mistral-embed")

```

```
pip install -qU langchain-nomic

```

```
import getpass
import os

if not os.environ.get("NOMIC_API_KEY"):
 os.environ["NOMIC_API_KEY"] = getpass.getpass("Enter API key for Nomic: ")

from langchain_nomic import NomicEmbeddings

embeddings = NomicEmbeddings(model="nomic-embed-text-v1.5")

```

```
pip install -qU langchain-nvidia-ai-endpoints

```

```
import getpass
import os

if not os.environ.get("NVIDIA_API_KEY"):
 os.environ["NVIDIA_API_KEY"] = getpass.getpass("Enter API key for NVIDIA: ")

from langchain_nvidia_ai_endpoints import NVIDIAEmbeddings

embeddings = NVIDIAEmbeddings(model="NV-Embed-QA")

```

```
pip install -qU langchain-voyageai

```

```
import getpass
import os

if not os.environ.get("VOYAGE_API_KEY"):
 os.environ["VOYAGE_API_KEY"] = getpass.getpass("Enter API key for Voyage AI: ")

from langchain-voyageai import VoyageAIEmbeddings

embeddings = VoyageAIEmbeddings(model="voyage-3")

```

```
pip install -qU langchain-ibm

```

```
import getpass
import os

if not os.environ.get("WATSONX_APIKEY"):
 os.environ["WATSONX_APIKEY"] = getpass.getpass("Enter API key for IBM watsonx: ")

from langchain_ibm import WatsonxEmbeddings

embeddings = WatsonxEmbeddings(
 model_id="ibm/slate-125m-english-rtrvr",
 url="https://us-south.ml.cloud.ibm.com",
 project_id="<WATSONX PROJECT_ID>",
)

```

```
pip install -qU langchain-core

```

```
from langchain_core.embeddings import DeterministicFakeEmbedding

embeddings = DeterministicFakeEmbedding(size=4096)

```

```
pip install -qU langchain-isaacus

```

```
import getpass
import os

if not os.environ.get("ISAACUS_API_KEY"):
os.environ["ISAACUS_API_KEY"] = getpass.getpass("Enter API key for Isaacus: ")

from langchain_isaacus import IsaacusEmbeddings

embeddings = IsaacusEmbeddings(model="kanon-2-embedder")

```

```
vector_1 = embeddings.embed_query(documents[0].page_content)
vector_2 = embeddings.embed_query(documents[1].page_content)

assert len(vector_1) == len(vector_2)
print(f"Generated vectors of length {len(vector_1)}\n")
print(vector_1[:10])

```

```
Generated vectors of length 1536

[-0.008586574345827103, -0.03341241180896759, -0.008936782367527485, -0.0036674530711025, 0.010564599186182022, 0.009598285891115665, -0.028587326407432556, -0.015824200585484505, 0.0030416189692914486, -0.012899317778646946]

```

Armed with a model for generating text embeddings, we can next store them in a special data structure that supports efficient similarity search.

## [​](#3-vector-stores)3. Vector stores

LangChain [VectorStore](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore) objects contain methods for adding text and [`Document`](https://reference.langchain.com/python/langchain-core/documents/base/Document) objects to the store, and querying them using various similarity metrics. They are often initialized with [embedding](/oss/python/integrations/embeddings) models, which determine how text data is translated to numeric vectors.
LangChain includes a suite of [integrations](/oss/python/integrations/vectorstores) with different vector store technologies. Some vector stores are hosted by a provider (e.g., various cloud providers) and require specific credentials to use; some (such as [Postgres](/oss/python/integrations/vectorstores/pgvector)) run in separate infrastructure that can be run locally or via a third-party; others can run in-memory for lightweight workloads. Let’s select a vector store:

{' ' * (self.list_depth - 1)}- In-memory
{' ' * (self.list_depth - 1)}- Amazon OpenSearch
{' ' * (self.list_depth - 1)}- AstraDB
{' ' * (self.list_depth - 1)}- Chroma
{' ' * (self.list_depth - 1)}- Milvus
{' ' * (self.list_depth - 1)}- MongoDB
{' ' * (self.list_depth - 1)}- PGVector
{' ' * (self.list_depth - 1)}- PGVectorStore
{' ' * (self.list_depth - 1)}- Pinecone
{' ' * (self.list_depth - 1)}- Qdrant
```
pip install -U "langchain-core"

```

```
from langchain_core.vectorstores import InMemoryVectorStore

vector_store = InMemoryVectorStore(embeddings)

```

```
pip install -qU boto3

```

```
from opensearchpy import RequestsHttpConnection

service = "es" # must set the service as 'es'
region = "us-east-2"
credentials = boto3.Session(
 aws_access_key_id="xxxxxx", aws_secret_access_key="xxxxx"
).get_credentials()
awsauth = AWS4Auth("xxxxx", "xxxxxx", region, service, session_token=credentials.token)

vector_store = OpenSearchVectorSearch.from_documents(
 docs,
 embeddings,
 opensearch_url="host url",
 http_auth=awsauth,
 timeout=300,
 use_ssl=True,
 verify_certs=True,
 connection_class=RequestsHttpConnection,
 index_name="test-index",
)

```

```
pip install -U "langchain-astradb"

```

```
from langchain_astradb import AstraDBVectorStore

vector_store = AstraDBVectorStore(
 embedding=embeddings,
 api_endpoint=ASTRA_DB_API_ENDPOINT,
 collection_name="astra_vector_langchain",
 token=ASTRA_DB_APPLICATION_TOKEN,
 namespace=ASTRA_DB_NAMESPACE,
)

```

```
pip install -qU langchain-chroma

```

```
from langchain_chroma import Chroma

vector_store = Chroma(
 collection_name="example_collection",
 embedding_function=embeddings,
 persist_directory="./chroma_langchain_db", # Where to save data locally, remove if not necessary
)

```

```
pip install -qU langchain-milvus

```

```
from langchain_milvus import Milvus

URI = "./milvus_example.db"

vector_store = Milvus(
 embedding_function=embeddings,
 connection_args={"uri": URI},
 index_params={"index_type": "FLAT", "metric_type": "L2"},
)

```

```
pip install -qU langchain-mongodb

```

```
from langchain_mongodb import MongoDBAtlasVectorSearch

vector_store = MongoDBAtlasVectorSearch(
 embedding=embeddings,
 collection=MONGODB_COLLECTION,
 index_name=ATLAS_VECTOR_SEARCH_INDEX_NAME,
 relevance_score_fn="cosine",
)

```

```
pip install -qU langchain-postgres

```

```
from langchain_postgres import PGVector

vector_store = PGVector(
 embeddings=embeddings,
 collection_name="my_docs",
 connection="postgresql+psycopg://...",
)

```

```
pip install -qU langchain-postgres

```

```
from langchain_postgres import PGEngine, PGVectorStore

pg_engine = PGEngine.from_connection_string(
 url="postgresql+psycopg://..."
)

vector_store = PGVectorStore.create_sync(
 engine=pg_engine,
 table_name='test_table',
 embedding_service=embeddings
)

```

```
pip install -qU langchain-pinecone

```

```
from langchain_pinecone import PineconeVectorStore
from pinecone import Pinecone

pc = Pinecone(api_key=...)
index = pc.Index(index_name)

vector_store = PineconeVectorStore(embedding=embeddings, index=index)

```

```
pip install -qU langchain-qdrant

```

```
from qdrant_client.models import Distance, VectorParams
from langchain_qdrant import QdrantVectorStore
from qdrant_client import QdrantClient

client = QdrantClient(":memory:")

vector_size = len(embeddings.embed_query("sample text"))

if not client.collection_exists("test"):
 client.create_collection(
 collection_name="test",
 vectors_config=VectorParams(size=vector_size, distance=Distance.COSINE)
 )
vector_store = QdrantVectorStore(
 client=client,
 collection_name="test",
 embedding=embeddings,
)

```

### [​](#seeding-the-vector-store)Seeding the vector store

Let’s seed the store with content from a PDF. [Here is a sample PDF](https://github.com/langchain-ai/langchain/blob/v0.3/docs/docs/example_data/nke-10k-2023.pdf) — a 10-k filing for Nike from 2023. We’ll read the PDF directly with a small helper and split it into smaller chunks before indexing.

```
import pypdf
from langchain_core.documents import Document

# Below is a minimal helper for demonstration purposes.
def load_pdf_pages(file_path: str) -> list[Document]:
 reader = pypdf.PdfReader(file_path)
 return [
 Document(
 page_content=page.extract_text() or "",
 metadata={"source": file_path, "page": i},
 )
 for i, page in enumerate(reader.pages)
 ]

file_path = "../example_data/nke-10k-2023.pdf"
docs = load_pdf_pages(file_path)
print(len(docs))

```

```
107

```

A page may be too coarse a representation for retrieval and downstream question-answering. Further splitting helps ensure that the meanings of relevant portions of the document are not “washed out” by surrounding text. We use [`RecursiveCharacterTextSplitter`](/oss/python/integrations/splitters), which recursively splits a document using common separators like new lines until each chunk is the appropriate size. This is the recommended text splitter for generic text use cases.
We set `add_start_index=True` so that the character index where each split Document starts within the initial Document is preserved as metadata attribute `start_index`.

```
from langchain_text_splitters import RecursiveCharacterTextSplitter

text_splitter = RecursiveCharacterTextSplitter(
 chunk_size=1000, chunk_overlap=200, add_start_index=True
)
all_splits = text_splitter.split_documents(docs)

print(len(all_splits))

```

```
514

```

We can now index the chunks into the vector store.

```
ids = vector_store.add_documents(documents=all_splits)

```

Note that most vector store implementations will allow you to connect to an existing vector store— e.g., by providing a client, index name, or other information. See the documentation for a specific [integration](/oss/python/integrations/vectorstores) for more detail.
Once we’ve instantiated a [`VectorStore`](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore) that contains documents, we can query it. [VectorStore](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore) includes methods for querying:

{' ' * (self.list_depth - 1)}- Synchronously and asynchronously;

{' ' * (self.list_depth - 1)}- By string query and by vector;

{' ' * (self.list_depth - 1)}- With and without returning similarity scores;

{' ' * (self.list_depth - 1)}- By similarity and [maximum marginal relevance](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore/max_marginal_relevance_search) (to balance similarity with query to diversity in retrieved results).

The methods will generally include a list of [Document](https://reference.langchain.com/python/langchain-core/documents/base/Document) objects in their outputs.
**Usage**
Embeddings typically represent text as a “dense” vector such that texts with similar meanings are geometrically close. This lets us retrieve relevant information just by passing in a question, without knowledge of any specific key-terms used in the document.
Return documents based on similarity to a string query:

```
results = vector_store.similarity_search(
 "How many distribution centers does Nike have in the US?"
)

print(results[0])

```

```
page_content='direct to consumer operations sell products through the following number of retail stores in the United States:
U.S. RETAIL STORES NUMBER
NIKE Brand factory stores 213
NIKE Brand in-line stores (including employee-only stores) 74
Converse stores (including factory stores) 82
TOTAL 369
In the United States, NIKE has eight significant distribution centers. Refer to Item 2. Properties for further information.
2023 FORM 10-K 2' metadata={'page': 4, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 3125}

```

Async query:

```
results = await vector_store.asimilarity_search("When was Nike incorporated?")

print(results[0])

```

```
page_content='Table of Contents
PART I
ITEM 1. BUSINESS
GENERAL
NIKE, Inc. was incorporated in 1967 under the laws of the State of Oregon. As used in this Annual Report on Form 10-K (this "Annual Report"), the terms "we," "us," "our,"
"NIKE" and the "Company" refer to NIKE, Inc. and its predecessors, subsidiaries and affiliates, collectively, unless the context indicates otherwise.
Our principal business activity is the design, development and worldwide marketing and selling of athletic footwear, apparel, equipment, accessories and services. NIKE is
the largest seller of athletic footwear and apparel in the world. We sell our products through NIKE Direct operations, which are comprised of both NIKE-owned retail stores
and sales through our digital platforms (also referred to as "NIKE Brand Digital"), to retail accounts and to a mix of independent distributors, licensees and sales' metadata={'page': 3, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 0}

```

Return scores:

```
# Note that providers implement different scores; the score here
# is a distance metric that varies inversely with similarity.

results = vector_store.similarity_search_with_score("What was Nike's revenue in 2023?")
doc, score = results[0]
print(f"Score: {score}\n")
print(doc)

```

```
Score: 0.23699893057346344

page_content='Table of Contents
FISCAL 2023 NIKE BRAND REVENUE HIGHLIGHTS
The following tables present NIKE Brand revenues disaggregated by reportable operating segment, distribution channel and major product line:
FISCAL 2023 COMPARED TO FISCAL 2022
•NIKE, Inc. Revenues were $51.2 billion in fiscal 2023, which increased 10% and 16% compared to fiscal 2022 on a reported and currency-neutral basis, respectively.
The increase was due to higher revenues in North America, Europe, Middle East & Africa ("EMEA"), APLA and Greater China, which contributed approximately 7, 6,
2 and 1 percentage points to NIKE, Inc. Revenues, respectively.
•NIKE Brand revenues, which represented over 90% of NIKE, Inc. Revenues, increased 10% and 16% on a reported and currency-neutral basis, respectively. This
increase was primarily due to higher revenues in Men's, the Jordan Brand, Women's and Kids' which grew 17%, 35%,11% and 10%, respectively, on a wholesale
equivalent basis.' metadata={'page': 35, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 0}

```

Return documents based on similarity to an embedded query:

```
embedding = embeddings.embed_query("How were Nike's margins impacted in 2023?")

results = vector_store.similarity_search_by_vector(embedding)
print(results[0])

```

```
page_content='Table of Contents
GROSS MARGIN
FISCAL 2023 COMPARED TO FISCAL 2022
For fiscal 2023, our consolidated gross profit increased 4% to $22,292 million compared to $21,479 million for fiscal 2022. Gross margin decreased 250 basis points to
43.5% for fiscal 2023 compared to 46.0% for fiscal 2022 due to the following:
*Wholesale equivalent
The decrease in gross margin for fiscal 2023 was primarily due to:
•Higher NIKE Brand product costs, on a wholesale equivalent basis, primarily due to higher input costs and elevated inbound freight and logistics costs as well as
product mix;
•Lower margin in our NIKE Direct business, driven by higher promotional activity to liquidate inventory in the current period compared to lower promotional activity in
the prior period resulting from lower available inventory supply;
•Unfavorable changes in net foreign currency exchange rates, including hedges; and
•Lower off-price margin, on a wholesale equivalent basis.
This was partially offset by:' metadata={'page': 36, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 0}

```

Learn more:

{' ' * (self.list_depth - 1)}- [API Reference](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore)

{' ' * (self.list_depth - 1)}- [Integration-specific docs](/oss/python/integrations/vectorstores)

## [​](#4-retrievers)4. Retrievers

LangChain [`VectorStore`](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStore) objects do not subclass [Runnable](https://reference.langchain.com/python/langchain-core/runnables/base/Runnable). LangChain [Retrievers](https://reference.langchain.com/python/langchain-core/retrievers/BaseRetriever) are Runnables, so they implement a standard set of methods (e.g., synchronous and asynchronous `invoke` and `batch` operations). Although we can construct retrievers from vector stores, retrievers can interface with non-vector store sources of data, as well (such as external APIs).
We can create a simple version of this ourselves, without subclassing `Retriever`. If we choose what method we wish to use to retrieve documents, we can create a runnable easily. Below we will build one around the `similarity_search` method:

```
from typing import List

from langchain_core.documents import Document
from langchain_core.runnables import chain

@chain
def retriever(query: str) -> List[Document]:
 return vector_store.similarity_search(query, k=1)

retriever.batch(
 [
 "How many distribution centers does Nike have in the US?",
 "When was Nike incorporated?",
 ],
)

```

```
[[Document(metadata={'page': 4, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 3125}, page_content='direct to consumer operations sell products through the following number of retail stores in the United States:\nU.S. RETAIL STORES NUMBER\nNIKE Brand factory stores 213 \nNIKE Brand in-line stores (including employee-only stores) 74 \nConverse stores (including factory stores) 82 \nTOTAL 369 \nIn the United States, NIKE has eight significant distribution centers. Refer to Item 2. Properties for further information.\n2023 FORM 10-K 2')],
 [Document(metadata={'page': 3, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 0}, page_content='Table of Contents\nPART I\nITEM 1. BUSINESS\nGENERAL\nNIKE, Inc. was incorporated in 1967 under the laws of the State of Oregon. As used in this Annual Report on Form 10-K (this "Annual Report"), the terms "we," "us," "our,"\n"NIKE" and the "Company" refer to NIKE, Inc. and its predecessors, subsidiaries and affiliates, collectively, unless the context indicates otherwise.\nOur principal business activity is the design, development and worldwide marketing and selling of athletic footwear, apparel, equipment, accessories and services. NIKE is\nthe largest seller of athletic footwear and apparel in the world. We sell our products through NIKE Direct operations, which are comprised of both NIKE-owned retail stores\nand sales through our digital platforms (also referred to as "NIKE Brand Digital"), to retail accounts and to a mix of independent distributors, licensees and sales')]]

```

Vectorstores implement an `as_retriever` method that will generate a Retriever, specifically a [`VectorStoreRetriever`](https://reference.langchain.com/python/langchain-core/vectorstores/base/VectorStoreRetriever). These retrievers include specific `search_type` and `search_kwargs` attributes that identify what methods of the underlying vector store to call, and how to parameterize them. For instance, we can replicate the above with the following:

```
retriever = vector_store.as_retriever(
 search_type="similarity",
 search_kwargs={"k": 1},
)

retriever.batch(
 [
 "How many distribution centers does Nike have in the US?",
 "When was Nike incorporated?",
 ],
)

```

```
[[Document(metadata={'page': 4, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 3125}, page_content='direct to consumer operations sell products through the following number of retail stores in the United States:\nU.S. RETAIL STORES NUMBER\nNIKE Brand factory stores 213 \nNIKE Brand in-line stores (including employee-only stores) 74 \nConverse stores (including factory stores) 82 \nTOTAL 369 \nIn the United States, NIKE has eight significant distribution centers. Refer to Item 2. Properties for further information.\n2023 FORM 10-K 2')],
 [Document(metadata={'page': 3, 'source': '../example_data/nke-10k-2023.pdf', 'start_index': 0}, page_content='Table of Contents\nPART I\nITEM 1. BUSINESS\nGENERAL\nNIKE, Inc. was incorporated in 1967 under the laws of the State of Oregon. As used in this Annual Report on Form 10-K (this "Annual Report"), the terms "we," "us," "our,"\n"NIKE" and the "Company" refer to NIKE, Inc. and its predecessors, subsidiaries and affiliates, collectively, unless the context indicates otherwise.\nOur principal business activity is the design, development and worldwide marketing and selling of athletic footwear, apparel, equipment, accessories and services. NIKE is\nthe largest seller of athletic footwear and apparel in the world. We sell our products through NIKE Direct operations, which are comprised of both NIKE-owned retail stores\nand sales through our digital platforms (also referred to as "NIKE Brand Digital"), to retail accounts and to a mix of independent distributors, licensees and sales')]]

```

`VectorStoreRetriever` supports search types of `"similarity"` (default), `"mmr"` (maximum marginal relevance, described above), and `"similarity_score_threshold"`. We can use the latter to threshold documents output by the retriever by similarity score.
Retrievers can easily be incorporated into more complex applications, such as [retrieval-augmented generation (RAG)](/oss/python/langchain/retrieval) applications that combine a given question with retrieved context into a prompt for a LLM. To learn more about building such an application, check out the [RAG tutorial](/oss/python/langchain/rag) tutorial.

## [​](#next-steps)Next steps

You’ve now seen how to build a semantic search engine over a PDF document.
For more on embeddings:

{' ' * (self.list_depth - 1)}- [Overview](/oss/python/langchain/retrieval)

{' ' * (self.list_depth - 1)}- [Available integrations](/oss/python/integrations/embeddings)

For more on vector stores:

{' ' * (self.list_depth - 1)}- [Overview](/oss/python/langchain/retrieval)

{' ' * (self.list_depth - 1)}- [Available integrations](/oss/python/integrations/vectorstores)

For more on RAG, see:

{' ' * (self.list_depth - 1)}- [Build a Retrieval Augmented Generation (RAG) App](/oss/python/langchain/rag)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/knowledge-base.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a data analysis agent from scratchPrevious](/oss/python/langchain/deep-agent-from-scratch)[Build a RAG agent with LangChainNext](/oss/python/langchain/rag)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
