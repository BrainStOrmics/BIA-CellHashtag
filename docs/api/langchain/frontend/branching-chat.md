# Branching chat

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Branching chat
Copy page

Edit messages, regenerate responses, and navigate conversation branchesCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Conversations with AI agents are rarely linear. You may want to rephrase a
question, regenerate a response you didn’t like, or explore a completely
different conversational path without losing your previous work. Branching chat
brings version-control semantics to your chat UI. Every edit creates a new
branch, and you can freely navigate between them.

This feature requires the [LangGraph Agent Server](/langsmith/local-server). Run your agent locally with `langgraph dev` or [deploy it to LangSmith](/langsmith/deployment) to use this pattern.

## [​](#what-is-branching-chat)What is branching chat?

Branching chat treats a conversation as a tree rather than a list. Each message
is a node, and editing a message or regenerating a response creates a **fork**
from that point. The original path is preserved as a sibling branch, so users
can switch back and forth between different conversation trajectories.
Key capabilities:

{' ' * (self.list_depth - 1)}- **Edit any user message:** rewrite a previous prompt and re-run the agent
from that point

{' ' * (self.list_depth - 1)}- **Regenerate any AI response:** ask the agent to produce a different answer
for the same input

{' ' * (self.list_depth - 1)}- **Navigate branches:** switch between different versions of the conversation
using per-message branch controls

## [​](#set-up-usestream-with-history)Set up useStream with history

To enable branching, pass `fetchStateHistory: true` so that `useStream`
retrieves checkpoint metadata needed for branch operations.
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
 assistantId: "branching_chat",
 fetchStateHistory: true,
 });

 return (
 <div>
 {stream.messages.map((msg) => {
 const metadata = stream.getMessagesMetadata(msg);
 return (
 <Message
 key={msg.id}
 message={msg}
 metadata={metadata}
 onEdit={(text) => handleEdit(stream, msg, metadata, text)}
 onRegenerate={() => handleRegenerate(stream, metadata)}
 onBranchSwitch={(id) => stream.setBranch(id)}
 />
 );
 })}
 </div>
 );
}

```

## [​](#understand-message-metadata)Understand message metadata

The `getMessagesMetadata(msg)` function returns branch information for each
message:

```
interface MessageMetadata {
 branch: string;
 branchOptions: string[];
 firstSeenState: {
 parent_checkpoint: Checkpoint | null;
 };
}

```

| 
 | Property | Description
 | `branch` | The branch ID of this specific message version
 | `branchOptions` | Array of all branch IDs available for this message position
 | `firstSeenState.parent_checkpoint` | The checkpoint just before this message. Use it as the fork point for edits and regenerations
When a message has only one version, `branchOptions` contains a single entry.
After an edit or regeneration, new branch IDs are added to `branchOptions`,
and you can navigate between them.

## [​](#edit-a-message)Edit a message

To edit a user message and create a new branch:

{' ' * (self.list_depth - 1)}- Get the `parent_checkpoint` from the message’s metadata

{' ' * (self.list_depth - 1)}- Submit the edited message with that checkpoint

{' ' * (self.list_depth - 1)}- The agent re-runs from that point, creating a new branch

```
function handleEdit(
 stream: ReturnType<typeof useStream>,
 originalMsg: HumanMessage,
 metadata: MessageMetadata,
 newText: string
) {
 const checkpoint = metadata.firstSeenState?.parent_checkpoint;
 if (!checkpoint) return;

 stream.submit(
 {
 messages: [{ ...originalMsg, content: newText }],
 },
 { checkpoint }
 );
}

```

After the edit:

{' ' * (self.list_depth - 1)}- The message’s `branchOptions` gains a new entry

{' ' * (self.list_depth - 1)}- The view switches to the new branch automatically

{' ' * (self.list_depth - 1)}- The agent re-runs from the fork point with the updated message

{' ' * (self.list_depth - 1)}- The original version is preserved and accessible via the branch switcher

## [​](#regenerate-a-response)Regenerate a response

To regenerate an AI response without changing the input:

{' ' * (self.list_depth - 1)}- Get the `parent_checkpoint` from the AI message’s metadata

{' ' * (self.list_depth - 1)}- Submit with `undefined` input and the parent checkpoint

{' ' * (self.list_depth - 1)}- The agent produces a fresh response, creating a new branch

```
function handleRegenerate(
 stream: ReturnType<typeof useStream>,
 metadata: MessageMetadata
) {
 const checkpoint = metadata.firstSeenState?.parent_checkpoint;
 if (!checkpoint) return;

 stream.submit(undefined, { checkpoint });
}

```

Each regeneration creates a new branch for the AI message at that position.
Users can then use the branch switcher to compare different responses.
Regeneration is useful for non-deterministic agents. Since LLM outputs vary
with temperature, regenerating the same prompt often produces meaningfully
different responses.

## [​](#build-a-branch-switcher)Build a branch switcher

When a message has multiple branches, show a compact inline control with the
current version index and navigation arrows:

```
function BranchSwitcher({
 metadata,
 onSwitch,
}: {
 metadata: MessageMetadata;
 onSwitch: (branchId: string) => void;
}) {
 const { branch, branchOptions } = metadata;

 if (branchOptions.length <= 1) return null;

 const currentIndex = branchOptions.indexOf(branch);
 const hasPrev = currentIndex > 0;
 const hasNext = currentIndex < branchOptions.length - 1;

 return (
 <div className="inline-flex items-center gap-1 rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-600">
 <button
 disabled={!hasPrev}
 onClick={() => onSwitch(branchOptions[currentIndex - 1])}
 className="hover:text-gray-900 disabled:opacity-30"
 aria-label="Previous version"
 >
 ◀
 </button>
 <span className="min-w-[3ch] text-center">
 {currentIndex + 1}/{branchOptions.length}
 </span>
 <button
 disabled={!hasNext}
 onClick={() => onSwitch(branchOptions[currentIndex + 1])}
 className="hover:text-gray-900 disabled:opacity-30"
 aria-label="Next version"
 >
 ▶
 </button>
 </div>
 );
}

```

When the user clicks a branch arrow, call `stream.setBranch(branchId)` to
switch the conversation view to that branch. This is instant since all branch
data is already loaded via `fetchStateHistory: true`.
Switching branches affects not only the target message but also all subsequent
messages. If you switch to a different version of message 3, messages 4, 5, 6,
etc. will also update to reflect the conversation that followed that version.

## [​](#how-branching-works-under-the-hood)How branching works under the hood

LangGraph persists every state transition as a **checkpoint**. When you submit
with a `checkpoint` parameter, the backend forks from that point instead of
appending to the current conversation. The result is a tree structure:

```
User: "What is React?"
 └─ AI: "React is a JavaScript library..." (branch A)
 └─ AI: "React is a UI framework..." (branch B, regenerated)

User: "Tell me about hooks" (branch A)
 └─ AI: "Hooks are functions..."

User: "Tell me about JSX" (edited from branch A)
 └─ AI: "JSX is a syntax extension..."

```

Each branch is an independent path through the conversation tree. Switching
branches updates the displayed messages but does not delete any data. All
branches persist in the checkpoint store.

## [​](#complete-message-component)Complete message component

Here is a full component that combines message display, editing, regeneration,
and branch switching:

```
function MessageWithBranching({
 message,
 metadata,
 stream,
}: {
 message: BaseMessage;
 metadata: MessageMetadata;
 stream: ReturnType<typeof useStream>;
}) {
 const [isEditing, setIsEditing] = useState(false);
 const [editText, setEditText] = useState(message.content as string);

 const isHuman = message._getType() === "human";
 const isAI = message._getType() === "ai";
 const hasBranches = metadata.branchOptions.length > 1;

 return (
 <div className="group relative py-2">
 {isEditing ? (
 <EditForm
 text={editText}
 onChange={setEditText}
 onSave={() => {
 handleEdit(stream, message as HumanMessage, metadata, editText);
 setIsEditing(false);
 }}
 onCancel={() => {
 setEditText(message.content as string);
 setIsEditing(false);
 }}
 />
 ) : (
 <>
 <div className={isHuman ? "text-right" : "text-left"}>
 <div
 className={
 isHuman
 ? "inline-block rounded-lg bg-blue-600 px-4 py-2 text-white"
 : "inline-block rounded-lg bg-gray-100 px-4 py-2"
 }
 >
 {message.content as string}
 </div>
 </div>

 <div className="mt-1 flex items-center gap-2 opacity-0 transition-opacity group-hover:opacity-100">
 {isHuman && (
 <button
 className="text-xs text-gray-400 hover:text-gray-700"
 onClick={() => setIsEditing(true)}
 >
 Edit
 </button>
 )}

 {isAI && (
 <button
 className="text-xs text-gray-400 hover:text-gray-700"
 onClick={() =>
 handleRegenerate(stream, metadata)
 }
 >
 Regenerate
 </button>
 )}

 {hasBranches && (
 <BranchSwitcher
 metadata={metadata}
 onSwitch={(id) => stream.setBranch(id)}
 />
 )}
 </div>
 </>
 )}
 </div>
 );
}

function EditForm({
 text,
 onChange,
 onSave,
 onCancel,
}: {
 text: string;
 onChange: (text: string) => void;
 onSave: () => void;
 onCancel: () => void;
}) {
 return (
 <div className="space-y-2">
 <textarea
 className="w-full rounded-lg border p-3 focus:outline-none focus:ring-2 focus:ring-blue-500"
 value={text}
 onChange={(e) => onChange(e.target.value)}
 rows={3}
 />
 <div className="flex gap-2">
 <button
 className="rounded bg-blue-600 px-4 py-1.5 text-sm text-white hover:bg-blue-700"
 onClick={onSave}
 >
 Save & Rerun
 </button>
 <button
 className="rounded border px-4 py-1.5 text-sm hover:bg-gray-50"
 onClick={onCancel}
 >
 Cancel
 </button>
 </div>
 </div>
 );
}

```

## [​](#combine-with-optimistic-updates)Combine with optimistic updates

Combine branching with optimistic updates for a seamless editing experience. When the user saves an edit, optimistically show the updated message before the server responds:

```
function handleEditOptimistic(
 stream: ReturnType<typeof useStream>,
 originalMsg: HumanMessage,
 metadata: MessageMetadata,
 newText: string
) {
 const checkpoint = metadata.firstSeenState?.parent_checkpoint;
 if (!checkpoint) return;

 const updatedMsg = { ...originalMsg, content: newText };

 stream.submit(
 { messages: [updatedMsg] },
 {
 checkpoint,
 optimisticValues: (prev) => {
 if (!prev?.messages) return { messages: [updatedMsg] };

 const idx = prev.messages.findIndex((m) => m.id === originalMsg.id);
 if (idx === -1) return prev;

 return {
 ...prev,
 messages: [...prev.messages.slice(0, idx), updatedMsg],
 };
 },
 }
 );
}

```

## [​](#add-keyboard-navigation)Add keyboard navigation

For power users, add keyboard shortcuts to navigate branches:

```
useEffect(() => {
 function handleKeyDown(e: KeyboardEvent) {
 if (!focusedMessageMetadata) return;

 const { branch, branchOptions } = focusedMessageMetadata;
 const idx = branchOptions.indexOf(branch);

 if (e.altKey && e.key === "ArrowLeft" && idx > 0) {
 stream.setBranch(branchOptions[idx - 1]);
 }
 if (e.altKey && e.key === "ArrowRight" && idx < branchOptions.length - 1) {
 stream.setBranch(branchOptions[idx + 1]);
 }
 }

 window.addEventListener("keydown", handleKeyDown);
 return () => window.removeEventListener("keydown", handleKeyDown);
}, [focusedMessageMetadata, stream]);

```

`Alt + ←` / `Alt + →` is a natural mapping for branch navigation since it
mirrors browser back/forward navigation.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Always enable `fetchStateHistory`**: without it, `getMessagesMetadata`
cannot return branch information.

{' ' * (self.list_depth - 1)}- **Only show the branch switcher when there are multiple branches**: a
`1/1` indicator adds clutter without value.

{' ' * (self.list_depth - 1)}- **Show branch controls on hover**: branch navigation arrows and edit buttons
should appear on hover to keep the UI clean.

{' ' * (self.list_depth - 1)}- **Keep the branch switcher compact**: it sits inline with message controls
and should not dominate the UI.

{' ' * (self.list_depth - 1)}- **Preserve scroll position**: when switching branches, try to keep the
viewport anchored to the message that changed.

{' ' * (self.list_depth - 1)}- **Indicate the active branch**: use subtle visual cues (e.g., a colored dot
or branch label) so users know which branch they’re viewing.

{' ' * (self.list_depth - 1)}- **Disable controls while streaming**: don’t allow edits or regeneration
while the agent is actively streaming a response. Check `stream.isLoading`
before enabling these actions.

{' ' * (self.list_depth - 1)}- **Preserve edit text on cancel**: if the user starts editing, then cancels,
reset the textarea to the original message content.

{' ' * (self.list_depth - 1)}- **Test with deep branch trees**: users who edit and regenerate frequently
can create many branches. Ensure the branch switcher and data handling
remain performant.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/branching-chat.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Human-in-the-LoopPrevious](/oss/python/langchain/frontend/human-in-the-loop)[Reasoning tokensNext](/oss/python/langchain/frontend/reasoning-tokens)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
