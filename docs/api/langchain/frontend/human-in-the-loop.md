# Human-in-the-Loop

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Human-in-the-Loop
Copy page

Add approval workflows with interrupt-based human reviewCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Not every agent action should run unsupervised. When an agent is about to send
an email, delete a record, execute a financial transaction, or perform any
irreversible operation, you need a human to review and approve the action first.
The Human-in-the-Loop (HITL) pattern lets your agent pause execution, present
the pending action to the user, and resume only after explicit approval.

## [​](#how-interrupts-work)How interrupts work

LangGraph agents support **interrupts**, explicit pause points where the agent
yields control back to the client. When the agent hits an interrupt:

{' ' * (self.list_depth - 1)}- The agent stops executing and emits an interrupt payload

{' ' * (self.list_depth - 1)}- The `useStream` hook surfaces the interrupt via `stream.interrupt`

{' ' * (self.list_depth - 1)}- Your UI renders a review card with approve/reject/edit options

{' ' * (self.list_depth - 1)}- The user makes a decision

{' ' * (self.list_depth - 1)}- Your code calls `stream.submit()` with a resume command

{' ' * (self.list_depth - 1)}- The agent picks up where it left off

## [​](#setting-up-usestream-for-hitl)Setting up useStream for HITL

Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";

const AGENT_URL = "http://localhost:2024";

export function Chat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "human_in_the_loop",
 });

 const interrupt = stream.interrupt;

 return (
 <div>
 {stream.messages.map((msg) => (
 <Message key={msg.id} message={msg} />
 ))}
 {interrupt && (
 <ApprovalCard
 interrupt={interrupt}
 onRespond={(response) =>
 stream.submit(null, { command: { resume: response } })
 }
 />
 )}
 </div>
 );
}

```

## [​](#the-interrupt-payload)The interrupt payload

When the agent pauses, `stream.interrupt` contains a `HITLRequest` with the
following structure:

```
interface HITLRequest {
 actionRequests: ActionRequest[];
 reviewConfigs: ReviewConfig[];
}

interface ActionRequest {
 action: string;
 args: Record<string, unknown>;
 description?: string;
}

interface ReviewConfig {
 allowedDecisions: ("approve" | "reject" | "edit" | "respond")[];
}

```

| 
 | Property | Description
 | `actionRequests` | Array of pending actions the agent wants to perform
 | `actionRequests[].action` | The action name (e.g. `"send_email"`, `"delete_record"`)
 | `actionRequests[].args` | Structured arguments for the action
 | `actionRequests[].description` | Optional human-readable description of what the action does
 | `reviewConfigs` | Per-action configuration controlling which decisions are allowed
 | `reviewConfigs[].allowedDecisions` | Which buttons to show: `"approve"`, `"reject"`, `"edit"`, `"respond"`

## [​](#decision-types)Decision types

The HITL pattern supports four decision types:

### [​](#approve)Approve

The user confirms the action should proceed as-is:

```
const response: HITLResponse = {
 decision: "approve",
};

stream.submit(null, { command: { resume: response } });

```

### [​](#reject)Reject

The user denies the action with an optional reason:

```
const response: HITLResponse = {
 decision: "reject",
 reason: "The email tone is too aggressive. Please revise.",
};

stream.submit(null, { command: { resume: response } });

```

When an action is rejected, the agent receives the rejection reason and can
decide how to proceed. It may rephrase, ask clarifying questions, or abandon
the action entirely.

### [​](#edit)Edit

The user modifies the action’s arguments before approving:

```
const response: HITLResponse = {
 decision: "edit",
 args: {
 ...originalArgs,
 subject: "Updated subject line",
 body: "Revised email body with softer language.",
 },
};

stream.submit(null, { command: { resume: response } });

```

### [​](#respond)Respond

The user provides a direct reply for “ask user” style tools. The `message` becomes the tool result and the tool itself is not executed:

```
const response: HITLResponse = {
 decision: "respond",
 message: "Blue.",
};

stream.submit(null, { command: { resume: response } });

```

Use `respond` when the tool is intentionally a placeholder for human input — for example, an `ask_user` tool that prompts the agent to collect information from the user.

## [​](#building-the-approvalcard)Building the ApprovalCard

Here is a full approval card component that handles all four decision types:

```
function ApprovalCard({
 interrupt,
 onRespond,
}: {
 interrupt: { value: HITLRequest };
 onRespond: (response: HITLResponse) => void;
}) {
 const request = interrupt.value;
 const [editedArgs, setEditedArgs] = useState(
 request.actionRequests[0]?.args ?? {}
 );
 const [rejectReason, setRejectReason] = useState("");
 const [respondMessage, setRespondMessage] = useState("");
 const [mode, setMode] = useState<"review" | "edit" | "reject" | "respond">("review");

 const action = request.actionRequests[0];
 const config = request.reviewConfigs[0];

 if (!action || !config) return null;

 return (
 <div className="rounded-lg border-2 border-amber-300 bg-amber-50 p-4">
 <h3 className="font-semibold text-amber-800">Action Review Required</h3>
 <p className="mt-1 text-sm text-amber-700">
 {action.description ?? `The agent wants to perform: ${action.action}`}
 </p>

 <div className="mt-3 rounded bg-white p-3 font-mono text-sm">
 <pre>{JSON.stringify(action.args, null, 2)}</pre>
 </div>

 {mode === "review" && (
 <div className="mt-4 flex gap-2">
 {config.allowedDecisions.includes("approve") && (
 <button
 className="rounded bg-green-600 px-4 py-2 text-white"
 onClick={() => onRespond({ decision: "approve" })}
 >
 Approve
 </button>
 )}
 {config.allowedDecisions.includes("reject") && (
 <button
 className="rounded bg-red-600 px-4 py-2 text-white"
 onClick={() => setMode("reject")}
 >
 Reject
 </button>
 )}
 {config.allowedDecisions.includes("edit") && (
 <button
 className="rounded bg-blue-600 px-4 py-2 text-white"
 onClick={() => setMode("edit")}
 >
 Edit
 </button>
 )}
 {config.allowedDecisions.includes("respond") && (
 <button
 className="rounded bg-purple-600 px-4 py-2 text-white"
 onClick={() => setMode("respond")}
 >
 Respond
 </button>
 )}
 </div>
 )}

 {mode === "reject" && (
 <div className="mt-4 space-y-2">
 <textarea
 className="w-full rounded border p-2"
 placeholder="Reason for rejection..."
 value={rejectReason}
 onChange={(e) => setRejectReason(e.target.value)}
 />
 <button
 className="rounded bg-red-600 px-4 py-2 text-white"
 onClick={() =>
 onRespond({ decision: "reject", reason: rejectReason })
 }
 >
 Confirm Rejection
 </button>
 </div>
 )}

 {mode === "edit" && (
 <div className="mt-4 space-y-2">
 <textarea
 className="w-full rounded border p-2 font-mono text-sm"
 value={JSON.stringify(editedArgs, null, 2)}
 onChange={(e) => {
 try {
 setEditedArgs(JSON.parse(e.target.value));
 } catch {
 // allow invalid JSON while editing
 }
 }}
 />
 <button
 className="rounded bg-blue-600 px-4 py-2 text-white"
 onClick={() =>
 onRespond({ decision: "edit", args: editedArgs })
 }
 >
 Submit Edits
 </button>
 </div>
 )}

 {mode === "respond" && (
 <div className="mt-4 space-y-2">
 <textarea
 className="w-full rounded border p-2"
 placeholder="Your response..."
 value={respondMessage}
 onChange={(e) => setRespondMessage(e.target.value)}
 />
 <button
 className="rounded bg-purple-600 px-4 py-2 text-white"
 onClick={() =>
 onRespond({ decision: "respond", message: respondMessage })
 }
 >
 Send Response
 </button>
 </div>
 )}
 </div>
 );
}

```

## [​](#the-resume-flow)The resume flow

After the user makes a decision, the full cycle looks like this:

{' ' * (self.list_depth - 1)}- Call `stream.submit(null, { command: { resume: hitlResponse } })`

{' ' * (self.list_depth - 1)}- The `useStream` hook sends the resume command to the LangGraph backend

{' ' * (self.list_depth - 1)}- The agent receives the `HITLResponse` and continues execution. The HITL response may be one of:

{' ' * (self.list_depth - 1)}- `"approve"`: The agent continues executing the next action

{' ' * (self.list_depth - 1)}- `"reject"`: The agent receives the rejection reasoning and decides its next step

{' ' * (self.list_depth - 1)}- `"edit"`: The agent runs the tool with the edited arguments

{' ' * (self.list_depth - 1)}- `"respond"`: The human’s message is returned directly as the tool result without executing the tool

{' ' * (self.list_depth - 1)}- The `interrupt` property resets to `null` as the agent resumes streaming

You can chain multiple HITL checkpoints in a single agent run. For example, an
agent might ask for approval to search, then ask again before sending an email
with the results. Each interrupt is handled independently.

## [​](#common-use-cases)Common use cases

| 
 | Use Case | Action | Review Config
 | Email sending | `send_email` | `["approve", "reject", "edit"]`
 | Database writes | `update_record` | `["approve", "reject"]`
 | Financial transactions | `transfer_funds` | `["approve", "reject"]`
 | File deletion | `delete_files` | `["approve", "reject"]`
 | API calls to external services | `call_api` | `["approve", "reject", "edit"]`
 | Collecting user input | `ask_user` | `["respond"]`

## [​](#handling-multiple-pending-actions)Handling multiple pending actions

An interrupt can contain multiple `actionRequests` when the agent wants to
perform several actions at once. Render a card for each and collect all
decisions before resuming:

```
function MultiActionReview({
 interrupt,
 onRespond,
}: {
 interrupt: { value: HITLRequest };
 onRespond: (responses: HITLResponse[]) => void;
}) {
 const [decisions, setDecisions] = useState<Record<number, HITLResponse>>({});
 const request = interrupt.value;

 const allDecided =
 Object.keys(decisions).length === request.actionRequests.length;

 return (
 <div className="space-y-4">
 {request.actionRequests.map((action, i) => (
 <SingleActionCard
 key={i}
 action={action}
 config={request.reviewConfigs[i]}
 onDecide={(response) =>
 setDecisions((prev) => ({ ...prev, [i]: response }))
 }
 />
 ))}
 {allDecided && (
 <button
 className="rounded bg-green-600 px-4 py-2 text-white"
 onClick={() =>
 onRespond(
 request.actionRequests.map((_, i) => decisions[i])
 )
 }
 >
 Submit All Decisions
 </button>
 )}
 </div>
 );
}

```

## [​](#best-practices)Best practices

Keep these guidelines in mind when implementing HITL workflows:

{' ' * (self.list_depth - 1)}- **Show clear context**. Always display *what* the agent wants to do and
*why*. Include the action description and the full arguments.

{' ' * (self.list_depth - 1)}- **Make approve the easiest path**. If the action looks correct, approving
should be a single click. Reserve multi-step flows for reject/edit.

{' ' * (self.list_depth - 1)}- **Validate edited args**. When users edit action arguments, validate the
JSON structure before sending. Show inline errors for malformed input.

{' ' * (self.list_depth - 1)}- **Persist the interrupt state**. If the user refreshes the page, the
interrupt should still be visible. `useStream` handles this via the thread’s
checkpoint.

{' ' * (self.list_depth - 1)}- **Log all decisions**. For audit trails, log every approve/reject/edit
decision with timestamps and the user who made the decision.

{' ' * (self.list_depth - 1)}- **Set timeouts thoughtfully**. Long-running agents should not block
indefinitely on human review. Consider showing how long the agent has been
waiting.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/human-in-the-loop.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Tool callingPrevious](/oss/python/langchain/frontend/tool-calling)[Branching chatNext](/oss/python/langchain/frontend/branching-chat)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
