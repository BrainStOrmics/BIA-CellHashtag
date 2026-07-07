# Human-in-the-loop

[Advanced usage](/oss/python/langchain/guardrails)

# Human-in-the-loop
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.The Human-in-the-Loop (HITL) [middleware](/oss/python/langchain/middleware/built-in#human-in-the-loop) lets you add human oversight to agent tool calls.
When a model proposes an action that might require review—for example, writing to a file or executing SQL—the middleware can pause execution and wait for a decision.
It does this by checking each tool call against a configurable policy. If intervention is needed, the middleware issues an [interrupt](https://reference.langchain.com/python/langgraph/types/interrupt) that halts execution. The graph state is saved using LangGraph’s [persistence layer](/oss/python/langgraph/persistence), so execution can pause safely and resume later.
A human decision then determines what happens next: the action can be approved as-is (`approve`), modified before running (`edit`), rejected with feedback (`reject`), or responded to directly (`respond`) for “ask user” style tools.

## [​](#interrupt-decision-types)Interrupt decision types

The [middleware](/oss/python/langchain/middleware/built-in#human-in-the-loop) defines four built-in ways a human can respond to an interrupt:

| 
 | Decision Type | Description | Example Use Case
 | ✅ `approve` | The action is approved as-is and executed without changes. | Send an email draft exactly as written
 | ✏️ `edit` | The tool call is executed with modifications. | Change the recipient before sending an email
 | ❌ `reject` | The tool call is rejected, with an explanation added to the conversation. | Reject an email draft and explain how to rewrite it
 | 💬 `respond` | Tool execution is skipped; the human’s message becomes the tool result. | Answer an “ask_user” prompt with a direct reply
The available decision types for each tool depend on the policy you configure in `interrupt_on`.
When multiple tool calls are paused at the same time, each action requires a separate decision.
Decisions must be provided in the same order as the actions appear in the interrupt request.
When **editing** tool arguments, make changes conservatively. Significant modifications to the original arguments may cause the model to re-evaluate its approach and potentially execute the tool multiple times or take unexpected actions.

## [​](#configuring-interrupts)Configuring interrupts

To use HITL, add the [middleware](/oss/python/langchain/middleware/built-in#human-in-the-loop) to the agent’s `middleware` list when creating the agent.
You configure it with a mapping of tool actions to the decision types that are allowed for each action. The middleware will interrupt execution when a tool call matches an action in the mapping.

```
from langchain.agents import create_agent
from langchain.agents.middleware import HumanInTheLoopMiddleware 
from langgraph.checkpoint.memory import InMemorySaver 

agent = create_agent(
 model="gpt-5.4",
 tools=[write_file, execute_sql, read_data],
 middleware=[
 HumanInTheLoopMiddleware(
 interrupt_on={
 "write_file": True, # All decisions (approve, edit, reject, respond) allowed
 "execute_sql": {"allowed_decisions": ["approve", "reject"]}, # No editing allowed
 "read_data": False, # Safe operation, no approval needed
 },
 # Prefix for interrupt messages - combined with tool name and args to form the full message
 # e.g., "Tool execution pending approval: execute_sql with query='DELETE FROM...'"
 # Individual tools can override this by specifying a "description" in their interrupt config
 description_prefix="Tool execution pending approval",
 ),
 ],
 # Human-in-the-loop requires checkpointing to handle interrupts.
 # In production, use a persistent checkpointer like AsyncPostgresSaver.
 checkpointer=InMemorySaver(),
)

```

You must configure a checkpointer to persist the graph state across interrupts.
In production, use a persistent checkpointer like [`AsyncPostgresSaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.postgres.aio.AsyncPostgresSaver). For testing or prototyping, use [`InMemorySaver`](https://reference.langchain.com/python/langgraph/checkpoints/#langgraph.checkpoint.memory.InMemorySaver).When invoking the agent, pass a `config` that includes the **thread ID** to associate execution with a conversation thread.
See the [LangGraph interrupts documentation](/oss/python/langgraph/interrupts) for details.

Configuration options[​](#param-interrupt-on)interrupt_ondictrequiredMapping of tool names to approval configs. Values can be `True` (interrupt with default config), `False` (auto-approve), or an `InterruptOnConfig` object.[​](#param-description-prefix)description_prefixstringdefault:"Tool execution requires approval"Prefix for action request descriptions**`InterruptOnConfig` options:**[​](#param-allowed-decisions)allowed_decisionslist[string]List of allowed decisions: `'approve'`, `'edit'`, `'reject'`, or `'respond'`[​](#param-description)descriptionstring | callableStatic string or callable function for custom description

## [​](#responding-to-interrupts)Responding to interrupts

When you invoke the agent, it runs until it either completes or an interrupt is raised. An interrupt is triggered when a tool call matches the policy you configured in `interrupt_on`. With `version="v2"`, the result is a `GraphOutput` with an `interrupts` attribute containing the actions that require review. You can then present those actions to a reviewer and resume execution once decisions are provided.

```
from langgraph.types import Command

# Human-in-the-loop leverages LangGraph's persistence layer.
# You must provide a thread ID to associate the execution with a conversation thread,
# so the conversation can be paused and resumed (as is needed for human review).
config = {"configurable": {"thread_id": "some_id"}}
# Run the graph until the interrupt is hit.
result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Delete old records from the database",
 }
 ]
 },
 config=config,
 version="v2",
)

# result is a GraphOutput with .value and .interrupts
print(result.interrupts)
# > (
# > Interrupt(
# > value={
# > 'action_requests': [
# > {
# > 'name': 'execute_sql',
# > 'arguments': {'query': 'DELETE FROM records WHERE created_at < NOW() - INTERVAL \'30 days\';'},
# > 'description': 'Tool execution pending approval\n\nTool: execute_sql\nArgs: {...}'
# > }
# > ],
# > 'review_configs': [
# > {
# > 'action_name': 'execute_sql',
# > 'allowed_decisions': ['approve', 'reject']
# > }
# > ]
# > }
# > ),
# > )

# Resume with approval decision
agent.invoke(
 Command(
 resume={"decisions": [{"type": "approve"}]} # or "reject"
 ),
 config=config, # Same thread ID to resume the paused conversation
 version="v2",
)

```

### [​](#decision-types)Decision types

{' ' * (self.list_depth - 1)}- ✅ approve
{' ' * (self.list_depth - 1)}- ✏️ edit
{' ' * (self.list_depth - 1)}- ❌ reject
{' ' * (self.list_depth - 1)}- 💬 respondUse `approve` to approve the tool call as-is and execute it without changes.
```
agent.invoke(
 Command(
 # Decisions are provided as a list, one per action under review.
 # The order of decisions must match the order of actions
 # in the interrupt request.
 resume={
 "decisions": [
 {
 "type": "approve",
 }
 ]
 }
 ),
 config=config, # Same thread ID to resume the paused conversation
 version="v2",
)

```
Use `edit` to modify the tool call before execution.
Provide the edited action with the new tool name and arguments.
```
agent.invoke(
 Command(
 # Decisions are provided as a list, one per action under review.
 # The order of decisions must match the order of actions
 # in the interrupt request.
 resume={
 "decisions": [
 {
 "type": "edit",
 # Edited action with tool name and args
 "edited_action": {
 # Tool name to call.
 # Will usually be the same as the original action.
 "name": "new_tool_name",
 # Arguments to pass to the tool.
 "args": {"key1": "new_value", "key2": "original_value"},
 }
 }
 ]
 }
 ),
 config=config, # Same thread ID to resume the paused conversation
 version="v2",
)

```
When **editing** tool arguments, make changes conservatively. Significant modifications to the original arguments may cause the model to re-evaluate its approach and potentially execute the tool multiple times or take unexpected actions.Use `reject` to reject the tool call and provide feedback instead of execution.
```
agent.invoke(
 Command(
 # Decisions are provided as a list, one per action under review.
 # The order of decisions must match the order of actions
 # in the interrupt request.
 resume={
 "decisions": [
 {
 "type": "reject",
 # An explanation about why the action was rejected
 "message": "No, this is wrong because ..., instead do this ...",
 }
 ]
 }
 ),
 config=config, # Same thread ID to resume the paused conversation
 version="v2",
)

```
The `message` is added to the conversation as feedback to help the agent understand why the action was rejected and what it should do instead.

---

### [​](#multiple-decisions)Multiple decisions
When multiple actions are under review, provide a decision for each action in the same order as they appear in the interrupt:
```
{
 "decisions": [
 {"type": "approve"},
 {
 "type": "edit",
 "edited_action": {
 "name": "tool_name",
 "args": {"param": "new_value"}
 }
 },
 {
 "type": "reject",
 "message": "This action is not allowed"
 }
 ]
}

```
Use `respond` for “ask user” style tools where the tool’s real implementation is the human’s reply. The `message` content is returned directly as the tool result; the tool itself is not executed.
```
agent.invoke(
 Command(
 # Decisions are provided as a list, one per action under review.
 # The order of decisions must match the order of actions
 # in the interrupt request.
 resume={
 "decisions": [
 {
 "type": "respond",
 # The human's reply, returned directly as the tool result
 "message": "Blue.",
 }
 ]
 }
 ),
 config=config, # Same thread ID to resume the paused conversation
 version="v2",
)

```
The `message` is returned to the agent as a successful `ToolMessage`. Use `respond` when the tool is intentionally a placeholder for human input—for example, an `ask_user` tool that prompts for clarification.

## [​](#streaming-with-human-in-the-loop)Streaming with human-in-the-loop

You can use `stream()` instead of `invoke()` to get real-time updates while the agent runs and handles interrupts. Use `stream_mode=['updates', 'messages']` with `version="v2"` to stream both agent progress and LLM tokens in the unified v2 format.

```
from langgraph.types import Command

config = {"configurable": {"thread_id": "some_id"}}

# Stream agent progress and LLM tokens until interrupt
for chunk in agent.stream(
 {"messages": [{"role": "user", "content": "Delete old records from the database"}]},
 config=config,
 stream_mode=["updates", "messages"],
 version="v2",
):
 if chunk["type"] == "messages":
 # LLM token
 token, metadata = chunk["data"]
 if token.content:
 print(token.content, end="", flush=True)
 elif chunk["type"] == "updates":
 # Check for interrupt
 if "__interrupt__" in chunk["data"]:
 print(f"\n\nInterrupt: {chunk['data']['__interrupt__']}")

# Resume with streaming after human decision
for chunk in agent.stream(
 Command(resume={"decisions": [{"type": "approve"}]}),
 config=config,
 stream_mode=["updates", "messages"],
 version="v2",
):
 if chunk["type"] == "messages":
 token, metadata = chunk["data"]
 if token.content:
 print(token.content, end="", flush=True)

```

See the [Streaming](/oss/python/langchain/streaming) guide for more details on stream modes.

## [​](#execution-lifecycle)Execution lifecycle

The middleware defines an `after_model` hook that runs after the model generates a response but before any tool calls are executed:

{' ' * (self.list_depth - 1)}- The agent invokes the model to generate a response.

{' ' * (self.list_depth - 1)}- The middleware inspects the response for tool calls.

{' ' * (self.list_depth - 1)}- If any calls require human input, the middleware builds a `HITLRequest` with `action_requests` and `review_configs` and calls [interrupt](https://reference.langchain.com/python/langgraph/types/interrupt).

{' ' * (self.list_depth - 1)}- The agent waits for human decisions.

{' ' * (self.list_depth - 1)}- Based on the `HITLResponse` decisions, the middleware executes approved or edited calls, synthesizes [ToolMessage](https://reference.langchain.com/python/langchain-core/messages/tool/ToolMessage)’s for rejected calls, returns human replies directly as [ToolMessage](https://reference.langchain.com/python/langchain-core/messages/tool/ToolMessage)’s for `respond` decisions, and resumes execution.

## [​](#custom-hitl-logic)Custom HITL logic

For more specialized workflows, you can build custom HITL logic directly using the [interrupt](https://reference.langchain.com/python/langgraph/types/interrupt) primitive and [middleware](/oss/python/langchain/middleware) abstraction.
Review the [execution lifecycle](#execution-lifecycle) above to understand how to integrate interrupts into the agent’s operation.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/human-in-the-loop.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Model Context Protocol (MCP)Previous](/oss/python/langchain/mcp)[Multi-agentNext](/oss/python/langchain/multi-agent)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
