# Skills

[Advanced usage](/oss/python/langchain/guardrails)[Multi-agent](/oss/python/langchain/multi-agent/index)

# Skills
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.In the **skills** architecture, specialized capabilities are packaged as invocable “skills” that augment an [agent’s](/oss/python/langchain/agents) behavior. Skills are primarily prompt-driven specializations that an agent can invoke on-demand.
For built-in skill support, see [Deep Agents](/oss/python/deepagents/skills).
This pattern is conceptually identical to [Agent Skills](https://agentskills.io/) and [llms.txt](https://llmstxt.org/) (introduced by Jeremy Howard), which uses tool calling for progressive disclosure of documentation. The skills pattern applies progressive disclosure to specialized prompts and domain knowledge rather than just documentation pages.For ready-to-use skills that improve your agent’s performance on LangChain ecosystem tasks, see the [LangChain Skills](https://github.com/langchain-ai/langchain-skills) repository.

## [​](#key-characteristics)Key characteristics

{' ' * (self.list_depth - 1)}- Prompt-driven specialization: Skills are primarily defined by specialized prompts

{' ' * (self.list_depth - 1)}- Progressive disclosure: Skills become available based on context or user needs

{' ' * (self.list_depth - 1)}- Team distribution: Different teams can develop and maintain skills independently

{' ' * (self.list_depth - 1)}- Lightweight composition: Skills are simpler than full sub-agents

{' ' * (self.list_depth - 1)}- Reference awareness: Skills can reference scripts, templates, and other resources

## [​](#when-to-use)When to use

Use the skills pattern when you want a single [agent](/oss/python/langchain/agents) with many possible specializations, you don’t need to enforce specific constraints between skills, or different teams need to develop capabilities independently. Common examples include coding assistants (skills for different languages or tasks), knowledge bases (skills for different domains), and creative assistants (skills for different formats).

## [​](#basic-implementation)Basic implementation

```
from langchain.tools import tool
from langchain.agents import create_agent

@tool
def load_skill(skill_name: str) -> str:
 """Load a specialized skill prompt.

 Available skills:
 - write_sql: SQL query writing expert
 - review_legal_doc: Legal document reviewer

 Returns the skill's prompt and context.
 """
 # Load skill content from file/database
 ...

agent = create_agent(
 model="gpt-5.4",
 tools=[load_skill],
 system_prompt=(
 "You are a helpful assistant. "
 "You have access to two skills: "
 "write_sql and review_legal_doc. "
 "Use load_skill to access them."
 ),
)

```

For a complete implementation, see the tutorial below.

## Tutorial: Build a SQL assistant with on-demand skills
Learn how to implement skills with progressive disclosure, where the agent loads specialized prompts and schemas on-demand rather than upfront.Learn more

## [​](#extending-the-pattern)Extending the pattern

When writing custom implementations, you can extend the basic skills pattern in several ways:

{' ' * (self.list_depth - 1)}- 
**Dynamic tool registration**: Combine progressive disclosure with state management to register new [tools](/oss/python/langchain/tools) as skills load. For example, loading a “database_admin” skill could both add specialized context and register database-specific tools (backup, restore, migrate). This uses the same tool-and-state mechanisms used across multi-agent patterns—tools updating state to dynamically change agent capabilities.

{' ' * (self.list_depth - 1)}- 
**Hierarchical skills**: Skills can define other skills in a tree structure, creating nested specializations. For instance, loading a “data_science” skill might make available sub-skills like “pandas_expert”, “visualization”, and “statistical_analysis”. Each sub-skill can be loaded independently as needed, allowing for fine-grained progressive disclosure of domain knowledge. This hierarchical approach helps manage large knowledge bases by organizing capabilities into logical groupings that can be discovered and loaded on-demand.

{' ' * (self.list_depth - 1)}- 
**Reference awareness**: While each skill only has one prompt, this prompt can reference the location of other assets and provide information on when the agent should use those assets.
When those assets become relevant, the agent will know that those files exist and read them into memory as needed to complete tasks.
This also follows the progressive disclosure pattern and limits the information in the context window.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/multi-agent/skills.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[HandoffsPrevious](/oss/python/langchain/multi-agent/handoffs)[RouterNext](/oss/python/langchain/multi-agent/router)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
